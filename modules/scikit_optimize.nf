#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process initialize_model {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/", mode: 'copy', overwrite: true
    cache 'lenient'
    fair true

    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed1), val(Rho_kg_per_m_cubed2), \
    path(statepoint1, stageAs: "statepoint1.json"), path(statepoint2, stageAs: "statepoint2.json"), \
    path(xsc1, stageAs: "xsc1.xsc"), path(xsc2, stageAs: "xsc2.xsc"),\
    path(coor1, stageAs: "coor1.coor"), path(coor2, stageAs: "coor2.coor"),\
    path(path_to_xml)
    output:
    tuple val(temp_K), path("initial_scikit_optimize_model.pkl"), val(0), emit: scikit_optimize_model
    path("log.txt")
    script:
    """
    #!/usr/bin/env python
    # for python3

    from typing import List
    from pydantic import BaseModel
    print("hello from ", $temp_K, $Rho_kg_per_m_cubed1,$Rho_kg_per_m_cubed2)

    class Point(BaseModel):
        density: float
        temperature: float
        pressure: float
        no_mol: float
        box_length: float
        rcut_couloumb: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.model_validate_json(json_data)

    loaded_point1 = load_point_from_json("$statepoint1")
    print("Loaded point1")
    print(loaded_point1)
    loaded_point2 = load_point_from_json("$statepoint2")
    print("Loaded point2")
    print(loaded_point2)

    liquid_box_length_Ang = loaded_point1.box_length
    vapor_box_length_Ang = loaded_point2.box_length
    percentage = 0.8
    RCC_START = 10.0
    RCC_END_BOX_0 = (float(liquid_box_length_Ang)/2.0)*percentage
    RCC_END_BOX_1 = (float(vapor_box_length_Ang)/2.0)*percentage

    import sys
    with open("log.txt", 'w') as sys.stdout:
        from skopt import Optimizer
        import pickle
        import numpy as np
        np.int = np.int64
        opt = Optimizer([(${params.alpha_lb},${params.alpha_ub}),\
                        (RCC_START,RCC_END_BOX_0),\
                        (${params.alpha_lb},${params.alpha_ub}),\
                        (RCC_START,RCC_END_BOX_1)],
                    "GP", acq_func="EI",
                    acq_optimizer="sampling",
                    initial_point_generator="lhs",
                    n_initial_points=${params.batch_size},
                    random_state=123)

        with open('initial_scikit_optimize_model.pkl', 'wb') as f:
            pickle.dump(opt, f)
        f.close()
        print("Model initialized with domain=",[(${params.alpha_lb},${params.alpha_ub}),\
                        (RCC_START,RCC_END_BOX_0),\
                        (${params.alpha_lb},${params.alpha_ub}),\
                        (RCC_START,RCC_END_BOX_1)])
        print("GP")
        print("acq_func=EI")
        print("acq_optimizer=sampling")
        print("initial_point_generator=lhs"),
        print("n_initial_points/batch size=",${params.batch_size})

    """
}


process build_two_box_system {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/ewald/input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed1), val(Rho_kg_per_m_cubed2), \
    path(statepoint1, stageAs: "statepoint1.json"), path(statepoint2, stageAs: "statepoint2.json"), \
    path(xsc1, stageAs: "xsc1.xsc"), path(xsc2, stageAs: "xsc2.xsc"),\
    path(coor1, stageAs: "coor1.coor"), path(coor2, stageAs: "coor2.coor"),\
    path(path_to_xml)

    output:
    tuple val(temp_K), path("system_liq.pdb"), path("system_liq.psf"), path("system_vap.pdb"), path("system_vap.psf"), path("system.inp"), \
    path(xsc1),path(coor1),path(xsc2),path(coor2), path("in_GEMC_NVT.conf"), emit: system
    script:
    """
    #!/usr/bin/env python
    from typing import List
    from pydantic import BaseModel
    print("hello from ", $temp_K, $Rho_kg_per_m_cubed1,$Rho_kg_per_m_cubed2)

    class Point(BaseModel):
        density: float
        temperature: float
        pressure: float
        no_mol: float
        box_length: float
        rcut_couloumb: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.model_validate_json(json_data)

    loaded_point1 = load_point_from_json("$statepoint1")
    print("Loaded point1")
    print(loaded_point1)
    loaded_point2 = load_point_from_json("$statepoint2")
    print("Loaded point2")
    print(loaded_point2)

    # GOMC Example for the NVT Ensemble using MoSDeF [1, 2, 5-10, 13-17]

    # Import the required packages and specify the force field (FF),
    # box dimensions, density, and mol ratios [1, 2, 5-10, 13-17].

    import mbuild as mb
    import unyt as u
    from mosdef_gomc.formats import charmm_writer as mf_charmm
    from mosdef_gomc.formats.charmm_writer import Charmm
    import mosdef_gomc.formats.gomc_conf_writer as gomc_control
    
    from   mbuild.lib.molecules.water import WaterSPC
    import foyer

    forcefield_file_water = "${path_to_xml}"

    liquid_box_length_Ang = loaded_point1.box_length

    liquid_box_density_kg_per_m_cubed = loaded_point1.density

    temperature = loaded_point1.temperature
    
    pressure=loaded_point1.pressure

    water = WaterSPC()
    water.name = 'SPCE'

    molecule_list = [water]
    residues_list = [water.name]
    fixed_bonds_angles_list = [water.name]

    ## Build the main liquid simulation box (box 0) for the simulation [1, 2, 13-17]

    liquid_water_box = mb.fill_box(compound=molecule_list,
                                        n_compounds=int(loaded_point1.no_mol),
                                        box=[liquid_box_length_Ang / 10,
                                            liquid_box_length_Ang / 10,
                                            liquid_box_length_Ang / 10]
                                        )

    vapor_box_length_Ang = loaded_point2.box_length

    vapor_box_density_kg_per_m_cubed = loaded_point2.density

    temperature = loaded_point1.temperature
    
    pressure=loaded_point1.pressure

    vapor_water_box = mb.fill_box(compound=molecule_list,
                                        n_compounds=int(loaded_point2.no_mol),
                                        box=[vapor_box_length_Ang / 10,
                                            vapor_box_length_Ang / 10,
                                            vapor_box_length_Ang / 10]
                                        )

    # Destroys water angles!!! Since GOMC is rigid water, dont use this.
    #water_box.energy_minimize(forcefield=forcefield_file_water , steps=10**5 )
    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    charmm = mf_charmm.Charmm(liquid_water_box,
                            'system_liq',
                            structure_box_1=vapor_water_box,
                            filename_box_1='system_vap',
                            ff_filename="system",
                            forcefield_selection=forcefield_file_water,
                            residues= residues_list,
                            gomc_fix_bonds_angles=fixed_bonds_angles_list
                            )


    import pickle
    # Unpickling the object
    with open("charmm.pkl", 'wb') as file:
         pickle.dump(charmm, file)

    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()

    from typing import List
    print("hello from ", $temp_K )

    from typing import Dict, Union
    from pydantic import BaseModel, Field
    import json

    class Point(BaseModel):
        density: float
        temperature: float
        pressure: float
        no_mol: float
        box_length: float
        rcut_couloumb: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.model_validate_json(json_data)

    loaded_point1 = load_point_from_json("$statepoint1")
    print("Loaded point1")
    print(loaded_point1)
    loaded_point2 = load_point_from_json("$statepoint2")
    print("Loaded point2")
    print(loaded_point2)

    import mbuild as mb
    import unyt as u
    from mosdef_gomc.formats import charmm_writer as mf_charmm
    from mosdef_gomc.formats.charmm_writer import Charmm
    import mosdef_gomc.formats.gomc_conf_writer as gomc_control

    # calc MC steps for gomc equilb
    # number of simulation steps
    if (${params.debugging}):
        gomc_steps_equilibration = 1000 #  set value for paper = 1 * 10**6
    else:
        gomc_steps_equilibration = 100000000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 10000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = int(gomc_steps_production/1000) # Monte Carlo Steps between console output
    coordinate_output_freq = int(gomc_steps_production/1000) # # set value for paper = 50 * 10**3
    restart_output_freq = int(gomc_steps_production/10) # # set value for paper = 50 * 10**3
    if (${params.debugging}):
        EqSteps = 100 # MCS for equilibration
        AdjSteps = 10 #MCS for adjusting max displacement, rotation, volume, etc.
    else:
        EqSteps = 100000 # MCS for equilibration
        AdjSteps = 1000 #MCS for adjusting max displacement, rotation, volume, etc.
    MC_steps = int(gomc_steps_production)
    # cutoff and tail correction
    Rcut_ang = 12 * u.angstrom
    Rcut_low_ang = 1.0 * u.angstrom
    LRC = True
    Exclude = "1-4"
    #RcutCoulomb_box_0=loaded_point1.rcut_couloumb
    #RcutCoulomb_box_1=loaded_point2.rcut_couloumb

    # MC move ratios
    DisFreq = 0.35
    RotFreq = 0.34
    VolFreq = 0.01
    MultiParticleFreq=0.00
    RegrowthFreq = 0.10
    IntraSwapFreq = 0.00
    IntraMEMC_2Freq = 0.00
    CrankShaftFreq = 0.00
    SwapFreq = 0.20
    MEMC_2Freq = 0.0

    gomc_output_data_every_X_steps = 50 * 10**3 # # set value for paper = 50 * 10**3
    # output all data and calc frequecy
    output_true_list_input = [
        True,
        int(gomc_output_data_every_X_steps),
    ]
    output_false_list_input = [
        False,
        int(gomc_output_data_every_X_steps),
    ]
    output_file_prefix="GOMC_GEMC_Production"


    gomc_control.write_gomc_control_file(charmm, 'in_GEMC_NVT.conf', 'GEMC_NVT', MC_steps, ${temp_K},
                                        Restart=True,
                                        check_input_files_exist=False,
                                        binCoordinates_box_0="${coor1}",
                                        extendedSystem_box_0="${xsc1}",
                                        binCoordinates_box_1="${coor2}",
                                        extendedSystem_box_1="${xsc2}",
                                        input_variables_dict={"VDWGeometricSigma": False,
                                                            "Ewald": True,
                                                            "ElectroStatic": True,
                                                            "PRNG": int(0),
                                                            "Pressure": None,
                                                            "Potential": "VDW",
                                                            "LRC": LRC,
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            #"RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            #"RcutCoulomb_box_1": RcutCoulomb_box_1,
                                                            "Exclude": Exclude,
                                                            "DisFreq": DisFreq,
                                                            "VolFreq": VolFreq,
                                                            "RotFreq": RotFreq,
                                                            "MultiParticleFreq": MultiParticleFreq,
                                                            "RegrowthFreq": RegrowthFreq,
                                                            "IntraSwapFreq": IntraSwapFreq,
                                                            "IntraMEMC-2Freq": IntraMEMC_2Freq,
                                                            "CrankShaftFreq": CrankShaftFreq,
                                                            "SwapFreq": SwapFreq,
                                                            "MEMC-2Freq": MEMC_2Freq,
                                                            "OutputName": output_file_prefix,
                                                            "EqSteps": EqSteps,
                                                            "AdjSteps":AdjSteps,
                                                            "PressureCalc": [True, pressure_calc_freq],
                                                            "RestartFreq": output_false_list_input,
                                                            "CheckpointFreq": output_false_list_input,
                                                            "DCDFreq": output_false_list_input,
                                                            "ConsoleFreq": [True, console_output_freq],
                                                            "BlockAverageFreq":[True, block_ave_output_freq],
                                                            "HistogramFreq": output_false_list_input,
                                                            "CoordinatesFreq": output_false_list_input,
                                                            "CBMC_First": 12,
                                                            "CBMC_Nth": 10,
                                                            "CBMC_Ang": 50,
                                                            "CBMC_Dih": 50,
                                                            }
                                        )
    """


}

process ask_points {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/batch/${iteration}/ask", mode: 'copy', overwrite: true
    debug false
    cache 'lenient'
    fair true
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf),\
    path(scikit_optimize_model), val(iteration)
    output:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf), val(iteration),\
    path("*.json"), emit: systems
    tuple val(temp_K), path("current_scikit_optimize_model.pkl"), val(iteration), emit: mdl
    tuple val(temp_K), path("*.json"), emit: json
    script:
    """
    #!/usr/bin/env python
    # for python3
    import sys
    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha_box_0: float
        rcc_box_0: float
        alpha_box_1: float
        rcc_box_1: float

    with open("log.txt", 'w') as sys.stdout:
        from skopt import Optimizer
        import pickle
        import numpy as np
        np.int = np.int64
        
        with open("${scikit_optimize_model}", 'rb') as f:
            opt = pickle.load(f)
        f.close()

        print("n_initial_points/batch size=",${params.batch_size})

        pointsList = opt.ask(n_points=int(${params.batch_size}))

        for count, value in enumerate(pointsList):
            # Create a Pydantic object
            point_obj = Point(alpha_box_0=value[0],\
                                rcc_box_0=value[1],\
                                alpha_box_1=value[2],\
                                rcc_box_1=value[3])

            # Serialize the Pydantic object to JSON
            with open(f"{count}.json", 'w') as file:
                file.write(point_obj.json())

        with open("current_scikit_optimize_model.pkl", 'wb') as f:
            pickle.dump(opt, f)
        f.close()

    """
}


process append_parameters_to_conf {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/batch/${iteration}/append_parameters_to_conf/${json.baseName}", mode: 'copy', overwrite: true
    cache 'lenient'
    fair true
    debug false
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf), val(iteration),\
    path(json)
    output:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf), val(iteration),\
    path(json), val(json.baseName), path("in_GEMC_NVT_RCC_ALPHA.conf")
    script:
    """
    #!/usr/bin/env python

    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha_box_0: float
        rcc_box_0: float
        alpha_box_1: float
        rcc_box_1: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.parse_raw(json_data)

    loaded_point = load_point_from_json("${json}")
    print("Loaded point")
    print(loaded_point)
    # Specify the file paths
    input_file_path = "${conf}"  # Replace with your actual input file path
    output_file_path = "in_GEMC_NVT_RCC_ALPHA.conf"  # Replace with your desired output file path
    print(input_file_path)
    print(output_file_path)
    # Open input file for reading and output file for writing
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        # Read the input file line by line
        for line in input_file:
            # Write each line to the output file
            output_file.write(line)
        print("WolfAlpha\t0\t{alpha}".format(alpha=loaded_point.alpha_box_0),file=output_file)
        print("RcutCoulomb\t0\t{alpha}".format(alpha=loaded_point.rcc_box_0),file=output_file)
        print("WolfAlpha\t1\t{alpha}".format(alpha=loaded_point.alpha_box_1),file=output_file)
        print("RcutCoulomb\t1\t{alpha}".format(alpha=loaded_point.rcc_box_1),file=output_file)

    """
}


process GOMC_GEMC_Production_Replica {
    // Important to allow for bad points
    // If GOMC doesn't finish, the objective function will be minimal.
    errorStrategy 'ignore'
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/batch/${iteration}/production/${json.baseName}", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf), val(iteration),\
    path(json), val(json_id), path(gemc_conf)
    output:
    tuple val(temp_K), val(iteration), val(json_id), path(json), path("GOMC_GEMC_Production.log"),  emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${gemc_conf} > local.conf
    GOMC_CPU_GEMC +p${task.cpus} local.conf > GOMC_GEMC_Production.log
    """
}


process Extract_Density_GOMC_GEMC_Production {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/batch/${iteration}/density/${json.baseName}", mode: 'copy', overwrite: true

    cpus 1
    debug false
    input:
    tuple val(temp_K), val(iteration), val(json_id), path(json), path(logfile)
    output:
    tuple val(temp_K),val(iteration), \
    path("${temp_K}_${iteration}_${json_id}_COMBINED.csv"),\
    emit: analysis
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import re
    import os
    import numpy as np
    from pathlib import Path
    def extract(filename,regex_pattern):

        EnRegex = re.compile(regex_pattern)
        print('extract_target',filename, 'pattern',regex_pattern)
        steps = []
        densities = []
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if EnRegex.match(line):
                        try:
                            steps.append(float(line.split()[1]))
                            densities.append(float(line.split()[8]))
                        except:
                            print(line)
                            print("An exception occurred") 
        except:
            print("Cant open",filename) 
        steps_np = np.array(steps)
        densities_np = np.array(densities)
        return steps_np, densities_np

    steps_box_0, densities_box_0 = extract("$logfile","STAT_0")
    df1 = pd.DataFrame(densities_box_0, index=None, columns=['${temp_K}_${iteration}_${json_id}_BOX_0'])
    #df.to_csv("${temp_K}_${iteration}_${json_id}_BOX_0.csv", header=True, sep=' ', index=False)

    steps_box_1, densities_box_1 = extract("$logfile","STAT_1")
    df2 = pd.DataFrame(densities_box_1, index=None, columns=['${temp_K}_${iteration}_${json_id}_BOX_1'])
    #df.to_csv("${temp_K}_${iteration}_${json_id}_BOX_1.csv", header=True, sep=' ', index=False)

    # Concatenate both dataframes
    df_combined = pd.concat([df1, df2], axis=1)

    # Save to CSV
    df_combined.to_csv("${temp_K}_${iteration}_${json_id}_COMBINED.csv", header=True, sep=' ', index=False)
    """
}


process Collate_GOMC_GEMC_Production { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/batch/${iteration}/merged/", mode: 'copy', overwrite: true
    cpus 1
    debug false
    input: tuple val(temp_K), val(iteration), path("*")
    output: path("merged_data.csv")

    shell:
    """
    paste * > merged_data.csv
    """
}


process tell_points {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize/temperature_${temp_K}_gemc/batch/${iteration}/ask", mode: 'copy', overwrite: true
    debug false
    cache 'lenient'
    fair true
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf),\
    path(scikit_optimize_model), val(iteration)
    output:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(conf), val(iteration),\
    path("*.json"), emit: systems
    tuple val(temp_K), path("current_scikit_optimize_model.pkl"), val(iteration), emit: mdl
    tuple val(temp_K), path("*.json"), emit: json
    script:
    """
    #!/usr/bin/env python
    # for python3
    import sys
    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha_box_0: float
        rcc_box_0: float
        alpha_box_1: float
        rcc_box_1: float

    with open("log.txt", 'w') as sys.stdout:
        from skopt import Optimizer
        import pickle
        import numpy as np
        np.int = np.int64
        
        with open("${scikit_optimize_model}", 'rb') as f:
            opt = pickle.load(f)
        f.close()

        print("n_initial_points/batch size=",${params.batch_size})

        pointsList = opt.ask(n_points=int(${params.batch_size}))

        for count, value in enumerate(pointsList):
            # Create a Pydantic object
            point_obj = Point(alpha_box_0=value[0],\
                                rcc_box_0=value[1],\
                                alpha_box_1=value[2],\
                                rcc_box_1=value[3])

            # Serialize the Pydantic object to JSON
            with open(f"{count}.json", 'w') as file:
                file.write(point_obj.json())

        with open("current_scikit_optimize_model.pkl", 'wb') as f:
            pickle.dump(opt, f)
        f.close()

    """
}

workflow initialize_scikit_optimize_model {
    take:
    system
    main:
    initialize_model(system)
    emit:
    scikit_optimize_model = initialize_model.out.scikit_optimize_model
    
}


workflow calibrate {
    take:
    scikit_optimize_model
    main:
    //ask_points(scikit_optimize_model)
    //ask_points.out.systems.transpose().view()
    //append_parameters_to_conf(ask_points.out.systems.transpose())

    //create_systems(ask.out.points.flatten())
    emit:
    scikit_optimize_model
    //scikit_optimize_model = ask_points.out.scikit_optimize_model
    //points = ask.out.points

}


workflow calibrate_wrapper {
    take:
    gemc_system_input
    ewald_density_data
    main:
    models = initialize_model(gemc_system_input)
    build_two_box_system(gemc_system_input)
    f = build_two_box_system.out.system.join(models.scikit_optimize_model)
    
    ask_points(f)
    append_parameters_to_conf(ask_points.out.systems.transpose())
    GOMC_GEMC_Production_Replica(append_parameters_to_conf.out)
    Extract_Density_GOMC_GEMC_Production(GOMC_GEMC_Production_Replica.out.record)
    collectedPoints = Extract_Density_GOMC_GEMC_Production.out.analysis.groupTuple(by:[0,1],size:params.batch_size,remainder:false,sort:true)
    sorted=collectedPoints.map { prefix1, prefix2, tbi -> tuple( prefix1, prefix2, tbi.sort{it.name}) }
    sorted.view()
    Collate_GOMC_GEMC_Production(sorted)
    //tellPointsInput = ask_points.out.mdl.join(collectedPoints, by:[0,1])
    //tellPointsInput.view()


}