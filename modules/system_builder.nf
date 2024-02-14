#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process build_solvent_system_from_torch {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${density}", mode: 'copy', overwrite: false

    debug false
    input:
    tuple val(density), path(statepoint), path(path_to_xml)
    output:
    tuple val(density), path("system.conf"), path("system.pdb"), path("system.psf"), path("system.inp"), emit: system

    script:
    """
    #!/usr/bin/env python

    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        density: float
        temperature: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.parse_raw(json_data)

    loaded_point = load_point_from_json("${statepoint}")
    print("Loaded point")
    print(loaded_point)
    print(loaded_point.density)
    print(loaded_point.temperature)
    density = loaded_point.density

    temperature = loaded_point.temperature

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

    liquid_box_length_Ang = 33

    liquid_box_density_kg_per_m_cubed = density

    water = WaterSPC()
    water.name = 'SPCE'

    molecule_list = [water]
    residues_list = [water.name]
    fixed_bonds_angles_list = [water.name]

    ## Build the main liquid simulation box (box 0) for the simulation [1, 2, 13-17]

    water_box = mb.fill_box(compound=molecule_list,
                                        density=liquid_box_density_kg_per_m_cubed,
                                        box=[liquid_box_length_Ang / 10,
                                            liquid_box_length_Ang / 10,
                                            liquid_box_length_Ang / 10]
                                        )
    # Destroys water angles!!! Since GOMC is rigid water, dont use this.
    #water_box.energy_minimize(forcefield=forcefield_file_water , steps=10**5 )
    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    charmm = mf_charmm.Charmm(water_box,
                            'system',
                            ff_filename="system",
                            forcefield_selection=forcefield_file_water,
                            residues= residues_list,
                            gomc_fix_bonds_angles=fixed_bonds_angles_list
                            )


    ## Write the write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    ### Note:  The electrostatics and Ewald are turned off in the
    # GOMC control file (i.e., False) since the n-alkanes beads in the
    # trappe-ua force field have no charge (i.e., the bead charges are all zero)
    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()

    #MC_steps=10000
    MC_steps=5000

    gomc_control.write_gomc_control_file(charmm, conf_filename='system',  ensemble_type='NVT', RunSteps=MC_steps, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": True,
                                                            "PRNG": int(0),
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1000],
                                                            "BlockAverageFreq":[False,1000],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "MultiParticleFreq":1.0,
                                                            "OutputName":"system"

                                                                }
                                        )




    """
}


process build_solvent_system {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/NVT/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), val(P_bar), val(No_mol), val(Rho_kg_per_m_cubed), val(L_m_if_cubed), val(RcutCoulomb), path(path_to_xml)
    output:
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path("statepoint.json"), path("system.pdb"), path("system.psf"), path("system.inp"), path("namd_system.inp"), emit: system
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path("statepoint.json"), emit: statepoint
    //path("system_npt.conf"), emit: npt_conf
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path("charmm.pkl"), emit: charmm

    script:
    """
    #!/usr/bin/env python

    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        density: float
        temperature: float
        pressure: float
        no_mol: float
        box_length: float
        rcut_couloumb: float

    # Create a Pydantic object
    point_obj = Point(density=${Rho_kg_per_m_cubed}, temperature=${temp_K}, pressure=${P_bar},no_mol=${No_mol},box_length=${L_m_if_cubed},rcut_couloumb=${RcutCoulomb})

    # Serialize the Pydantic object to JSON
    with open("statepoint.json", 'w') as file:
        file.write(point_obj.model_dump_json())

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

    liquid_box_length_Ang = ${L_m_if_cubed}

    liquid_box_density_kg_per_m_cubed = ${Rho_kg_per_m_cubed}

    temperature = ${temp_K}
    
    pressure=${P_bar}

    water = WaterSPC()
    water.name = 'SPCE'

    molecule_list = [water]
    residues_list = [water.name]
    fixed_bonds_angles_list = [water.name]

    ## Build the main liquid simulation box (box 0) for the simulation [1, 2, 13-17]

    water_box = mb.fill_box(compound=molecule_list,
                                        #density=${Rho_kg_per_m_cubed},
                                        n_compounds=int(${No_mol}),
                                        box=[liquid_box_length_Ang / 10,
                                            liquid_box_length_Ang / 10,
                                            liquid_box_length_Ang / 10]
                                        )
    # Destroys water angles!!! Since GOMC is rigid water, dont use this.
    #water_box.energy_minimize(forcefield=forcefield_file_water , steps=10**5 )
    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    charmm = mf_charmm.Charmm(water_box,
                            'system',
                            ff_filename="system",
                            forcefield_selection=forcefield_file_water,
                            residues= residues_list,
                            gomc_fix_bonds_angles=fixed_bonds_angles_list
                            )

    namd_charmm = mf_charmm.Charmm(water_box,
                            'namd_system',
                            ff_filename="namd_system",
                            forcefield_selection=forcefield_file_water,
                            residues= residues_list,
                            gomc_fix_bonds_angles=None
                            )

    ## Write the write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]
    namd_charmm.write_inp()

    ### Note:  The electrostatics and Ewald are turned off in the
    # GOMC control file (i.e., False) since the n-alkanes beads in the
    # trappe-ua force field have no charge (i.e., the bead charges are all zero)
    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()

    import pickle
    # Unpickling the object
    with open("charmm.pkl", 'wb') as file:
         pickle.dump(charmm, file)

    """
}


process build_two_box_system_already_calibrated {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/wolf/input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed1), val(Rho_kg_per_m_cubed2), \
    path(statepoint1, stageAs: "statepoint1.json"), path(statepoint2, stageAs: "statepoint2.json"), \
    path(xsc1, stageAs: "xsc1.xsc"), path(xsc2, stageAs: "xsc2.xsc"),\
    path(coor1, stageAs: "coor1.coor"), path(coor2, stageAs: "coor2.coor"),\
    path(pdb1, stageAs: "psb1.pdb"), path(pdb2, stageAs: "pdb2.pdb"),\
    path(psf1, stageAs: "psf1.psf"), path(psf2, stageAs: "psf2.psf"),\
    path(chk, stageAs: "chk.chk"), path(convergence_obj),\
    path(path_to_xml), val(METHOD)

    output:
    output:
    tuple val(temp_K), val(METHOD), path(pdb1), path(psf1), path(pdb2), path(psf2), path("system.inp"), \
    path(xsc1),path(coor1),path(xsc2),path(coor2),path(chk),path("in_GEMC_NVT.conf"), emit: system
    tuple val(temp_K), path("charmm.pkl"), path(xsc1),path(coor1),path(xsc2),path(coor2),path(statepoint1),path(statepoint2), emit: charmm
    tuple val(temp_K), path("charmm.pkl"),path(statepoint1),path(statepoint2), emit: charmm_norestarts
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


    # Convergence data
    from typing import Dict, Union
    from pydantic import BaseModel, Field
    import json

    class InnerModel(BaseModel):
        ConvergedRCut: float
        ConvergedAlpha: float
        ConvergedFAlpha: float
        ConvergedDFDAlpha: float

    class FooBarModel(BaseModel):
        models: Dict[str, Dict[int, InnerModel]]  # Adjusted for the nested dictionary

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> FooBarModel:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return FooBarModel.model_validate_json(json_data)
    
    convergence_obj = load_point_from_json("${convergence_obj}")

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

    # calc MC steps for gomc equilb
    # number of simulation steps
    if (${params.debugging}):
        gomc_steps_equilibration = 1000 #  set value for paper = 1 * 10**6
    else:
        gomc_steps_equilibration = 10000000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 10000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = 100000 # Monte Carlo Steps between console output
    coordinate_output_freq = 100000 # # set value for paper = 50 * 10**3
    restart_output_freq = 100000 # # set value for paper = 50 * 10**3
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
    RcutCoulomb_box_0=convergence_obj.models["${METHOD}"][0].ConvergedRCut
    RcutCoulomb_box_1=convergence_obj.models["${METHOD}"][1].ConvergedRCut

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
                                        Coordinates_box_0="${pdb1}",
                                        Coordinates_box_1="${pdb2}",
                                        Structure_box_0="${psf1}",
                                        Structure_box_1="${psf2}",
                                        binCoordinates_box_0="${coor1}",
                                        extendedSystem_box_0="${xsc1}",
                                        binCoordinates_box_1="${coor2}",
                                        extendedSystem_box_1="${xsc2}",
                                        input_variables_dict={"VDWGeometricSigma": False,
                                                            "Ewald": False,
                                                            "ElectroStatic": True,
                                                            "PRNG": int(0),
                                                            "Pressure": None,
                                                            "Potential": "VDW",
                                                            "LRC": LRC,
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            "RcutCoulomb_box_1": RcutCoulomb_box_1,
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
                                                            "RestartFreq": [True, restart_output_freq],
                                                            "CheckpointFreq": [True, restart_output_freq],
                                                            "DCDFreq": [True, coordinate_output_freq],
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

    kind_pot = "${METHOD}".split('_')
    file1 = open("in_GEMC_NVT.conf", "a")
    defAlphaLine = "{box}\\t{val}\\n".format(box="Wolf", val="True")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\n".format(box="WolfKind", val=kind_pot[0])
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\n".format(box="WolfPotential", val=kind_pot[1])
    file1.writelines(defAlphaLine)
    defAlphaLine = "{key}\\t{box}\\t{val}\\n".format(key="WolfAlpha",box="0", val=convergence_obj.models["${METHOD}"][0].ConvergedAlpha)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{key}\\t{box}\\t{val}\\n".format(key="WolfAlpha",box="1", val=convergence_obj.models["${METHOD}"][1].ConvergedAlpha)
    file1.writelines(defAlphaLine)

    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file="${chk}")
    file1.writelines(defAlphaLine)
    """
}


process build_two_box_system {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/ewald/input", mode: 'copy', overwrite: false
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
    path(xsc1),path(coor1),path(xsc2),path(coor2),path("charmm.pkl"),path(statepoint1),path(statepoint2), emit: system
    tuple val(temp_K), path("charmm.pkl"), path(xsc1),path(coor1),path(xsc2),path(coor2),path(statepoint1),path(statepoint2), emit: charmm
    tuple val(temp_K), path("charmm.pkl"),path(statepoint1),path(statepoint2), emit: charmm_norestarts
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

    """
}


process build_two_box_system_calibrate {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/calibration/input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed1), val(Rho_kg_per_m_cubed2), \
    path(statepoint1, stageAs: "statepoint1.json"), path(statepoint2, stageAs: "statepoint2.json"), \
    path(xsc1, stageAs: "xsc1.xsc"), path(xsc2, stageAs: "xsc2.xsc"),\
    path(coor1, stageAs: "coor1.coor"), path(coor2, stageAs: "coor2.coor"),\
    path(pdb1, stageAs: "psb1.pdb"), path(pdb2, stageAs: "pdb2.pdb"),\
    path(psf1, stageAs: "psf1.psf"), path(psf2, stageAs: "psf2.psf"),\
    path(chk, stageAs: "chk.chk"),\
    path(path_to_xml)

    output:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path("system.inp"), \
    path(xsc1),path(coor1),path(xsc2),path(coor2),path(chk),path("in_GEMC_NVT.conf"), emit: system
    tuple val(temp_K), path("charmm.pkl"), path(xsc1),path(coor1),path(xsc2),path(coor2),path(statepoint1),path(statepoint2), emit: charmm
    tuple val(temp_K), path("charmm.pkl"),path(statepoint1),path(statepoint2), emit: charmm_norestarts
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

    # calc MC steps for gomc equilb
    # number of simulation steps
    if (${params.debugging}):
        gomc_steps_equilibration = 1000 #  set value for paper = 1 * 10**6
    else:
        gomc_steps_equilibration = 100000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 10000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = 100000 # Monte Carlo Steps between console output
    coordinate_output_freq = 100000 # # set value for paper = 50 * 10**3
    restart_output_freq = 10000 # # set value for paper = 50 * 10**3
    if (${params.debugging}):
        EqSteps = 100 # MCS for equilibration
        AdjSteps = 10 #MCS for adjusting max displacement, rotation, volume, etc.
    else:
        EqSteps = 1000 # MCS for equilibration
        AdjSteps = 100 #MCS for adjusting max displacement, rotation, volume, etc.
    MC_steps = int(gomc_steps_production)
    # cutoff and tail correction
    Rcut_ang = 12 * u.angstrom
    Rcut_low_ang = 1.0 * u.angstrom
    LRC = True
    Exclude = "1-4"
    RcutCoulomb_box_0=loaded_point1.rcut_couloumb
    RcutCoulomb_box_1=loaded_point2.rcut_couloumb

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
                                        Coordinates_box_0="${pdb1}",
                                        Coordinates_box_1="${pdb2}",
                                        Structure_box_0="${psf1}",
                                        Structure_box_1="${psf2}",
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
                                                            "RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            "RcutCoulomb_box_1": RcutCoulomb_box_1,
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
                                                            "RestartFreq": [True, restart_output_freq],
                                                            "CheckpointFreq": [True, restart_output_freq],
                                                            "DCDFreq": [True, coordinate_output_freq],
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
    if (${params.debugging}):
        NUM_POINTS = 5.0
    else:
        NUM_POINTS = 50.0
    percentage = 0.80
    ALPHA_START = 0.0
    ALPHA_END = 0.5
    ALPHA_DELTA = (ALPHA_END-ALPHA_START)/(NUM_POINTS-1)
    RCC_START = 10.0
    RCC_END_BOX_0 = (float(liquid_box_length_Ang)/2.0)*percentage
    RCC_DELTA_BOX_0 = (RCC_END_BOX_0-RCC_START)/(NUM_POINTS-1)
    RCC_END_BOX_1 = (float(vapor_box_length_Ang)/2.0)*percentage
    RCC_DELTA_BOX_1 = (RCC_END_BOX_1-RCC_START)/(NUM_POINTS-1)
    file1 = open("in_GEMC_NVT.conf", "a")
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="10000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start=RCC_START,\
    end=RCC_END_BOX_0,delta=RCC_DELTA_BOX_0)
    file1.writelines(defAlphaLine)

    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="1",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="1",start=RCC_START,\
    end=RCC_END_BOX_1,delta=RCC_DELTA_BOX_1)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file="${chk}")
    file1.writelines(defAlphaLine)
    """
}


process write_namd_confs {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/NVT/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/namd_input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path("statepoint.json"), path("system.pdb"), path("system.psf"), path("system.inp"), path("namd_system.inp")
    tuple path(path_to_minimization_template), path(path_to_nvt_template), path(path_to_npt_template)
    output:
    //tuple val(temp_K), val(Rho_kg_per_m_cubed),path("statepoint.json"),path("system_nvt.conf"), path("system_npt.conf"), path("system.pdb"), path("system.psf"), path("system.inp"), emit: system
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path("namd_min.conf"), path("namd_nvt.conf"), path("namd_npt.conf"), emit: namd

    script:
    """
    #!/usr/bin/env python

    import os
    def _setup_conf(fname, template, data, overwrite=False):

        #Create conf files based on a template and provided data.

        #Parameters
        #----------
        #fname: str
        #    Name of the file to be saved out
        #template: str, or jinja2.Template
        #    Either a jinja2.Template or path to a jinja template
        #data: dict
        #    Dictionary storing data matched with the fields available in the template
        #overwrite: bool, optional, default=False
        #    Options to overwrite (or not) existing mdp file of the

        #Returns
        #-------
        #File saved with names defined by fname
        

        from jinja2 import Template

        if isinstance(template, str):
            with open(template, "r") as f:
                template = Template(f.read())

        if not overwrite:
            if os.path.isfile(fname):
                raise FileExistsError(
                    f"{fname} already exists. Set overwrite=True to write out."
                )

        rendered = template.render(data)
        with open(fname, "w") as f:
            f.write(rendered)

        return None

    from typing import List
    from pydantic import BaseModel

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

    loaded_point = load_point_from_json("statepoint.json")
    print("Loaded point")
    print(loaded_point)


    coords = {}

    #minsteps=10000
    #nvt_eq_steps=50000
    #npt_eq_steps=100000
    # gas phase needs less equilibration
    if (${Rho_kg_per_m_cubed}<400):
        if (${params.debugging}):
            minsteps=500
            nvt_eq_steps=500
            npt_eq_steps=500
        else:
            minsteps=500
            nvt_eq_steps=5000
            npt_eq_steps=5000
    else:
        if (${params.debugging}):
            minsteps=500
            nvt_eq_steps=1000
            npt_eq_steps=1000
        else:
            minsteps=500
            nvt_eq_steps=20000
            npt_eq_steps=20000

    coords["waterModel"]="TIP3"
    coords["temp"]=loaded_point.temperature
    coords["rcut_couloumb"]=loaded_point.rcut_couloumb
    coords["pairlistdist"]=loaded_point.rcut_couloumb+(loaded_point.rcut_couloumb*0.1)
    coords["outputname"]="minimization"
    coords["minimize_steps"]=minsteps

    coords["structure"]="system.psf"
    coords["coordinates"]="system.pdb"

    coords["gomcwaterparameters"]="system.inp"
    coords["namdwaterparameters"]="namd_system.inp"

    coords["X_DIM_BOX"]=loaded_point.box_length
    coords["Y_DIM_BOX"]=loaded_point.box_length
    coords["Z_DIM_BOX"]=loaded_point.box_length
    coords["alpha"]=90
    coords["beta"]=90
    coords["gamma"]=90
    coords["X_ORIGIN_BOX"]=loaded_point.box_length/2
    coords["Y_ORIGIN_BOX"]=loaded_point.box_length/2
    coords["Z_ORIGIN_BOX"]=loaded_point.box_length/2
    _setup_conf("namd_min.conf", "$path_to_minimization_template", coords, overwrite=False)
    coords["eq_steps"]=nvt_eq_steps
    coords["bincoordinates"]="min.restart.coor"
    coords["outputname"]="nvt_eq"
    _setup_conf("namd_nvt.conf", "$path_to_nvt_template", coords, overwrite=False)
    coords["eq_steps"]=npt_eq_steps
    coords["bincoordinates"]="nvt_equil.restart.coor"
    coords["pressure"]=loaded_point.pressure
    coords["outputname"]="npt_eq"
    _setup_conf("namd_npt.conf", "$path_to_npt_template", coords, overwrite=False)
    """
}


process NAMD_equilibration_solvent_system {
    cache 'lenient'
    fair true
    container "${params.container__namd}"
    publishDir "${params.output_folder}/systems/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/namd_npt_eq", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint),path(pdb), path(psf), path(inp), path(namd_inp)
    tuple path(namd_minimization_conf), path(namd_nvt_conf), path(namd_npt_conf)
    output:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint), path(pdb), path(psf), path(inp), path("npt_equil.restart.xsc"), path("npt_equil.restart.coor"), emit: system
    tuple path("npt_equil.restart.xsc"), path("npt_equil.restart.coor"), emit: restart_files
    tuple path("minimization.log"), path("nvt_equil.log"), path("npt_equil.log"),  emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${namd_minimization_conf} > local.conf
    namd2 +p${task.cpus} local.conf > minimization.log
    cat ${namd_nvt_conf} > local.conf
    namd2 +p${task.cpus} local.conf > nvt_equil.log
    cat ${namd_npt_conf} > local.conf
    namd2 +p${task.cpus} local.conf > npt_equil.log
    """
}


process NAMD_NVT_equilibration {
    cache 'lenient'
    fair true
    container "${params.container__namd}"
    publishDir "${params.output_folder}/NVT/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/namd_nvt_eq", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed),val(temp_K),path(statepoint),path(pdb), path(psf), path(inp), path(namd_inp),val(temp_K),path(namd_minimization_conf), path(namd_nvt_conf), path(namd_npt_conf)
    output:
    tuple val(temp_K), val(Rho_kg_per_m_cubed), path(statepoint), path("nvt_equil.restart.xsc"), path("nvt_equil.restart.coor"), emit: restart_files
    tuple path("minimization.log"), path("nvt_equil.log"), emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${namd_minimization_conf} > local.conf
    namd2 +p${task.cpus} local.conf > minimization.log
    cat ${namd_nvt_conf} > local.conf
    namd2 +p${task.cpus} local.conf > nvt_equil.log
    """
}


process GOMC_NPT_equilibration {
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_npt_eq", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint),path(pdb), path(psf), path(inp), path(namd_inp)
    path(npt_conf)
    tuple path(restart_xsc), path(restart_coor)
    output:
    tuple path("system_npt_BOX_0_restart.xsc"), path("system_npt_BOX_0_restart.coor"), path("system_npt_restart.chk"), emit: restart_files
    tuple path("NPT_equilibration.log"), path(npt_conf),  emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${npt_conf} > local.conf
    GOMC_CPU_NPT +p${task.cpus} local.conf > NPT_equilibration.log
    """
}


process NVT_equilibration_solvent_system {
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_nvt_eq", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint),path(nvt_conf), path(npt_conf), path(pdb), path(psf), path(inp)
    output:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint), path(npt_conf), path(pdb), path(psf), path(inp),\
    path("system_nvt_restart.chk"),path("system_nvt_BOX_0_restart.coor"), path("system_nvt_BOX_0_restart.xsc"), emit: system
    tuple path("NVT_equilibration.log"), path(nvt_conf),  emit: record

    shell:
    """
    
    #!/bin/bash
    cat ${nvt_conf} > local.conf
    GOMC_CPU_NVT +p${task.cpus} local.conf > NVT_equilibration.log
    """
}


process write_gomc_confs {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_npt_eq_input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(json), path(charmm)
    tuple path(restart_xsc), path(restart_coor)
    output:
    path("system_npt.conf"), emit: npt_conf

    script:
    """
    #!/usr/bin/env python
    from typing import List
    from pydantic import BaseModel

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

    loaded_point = load_point_from_json("${json}")

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

    liquid_box_length_Ang = loaded_point.box_length

    liquid_box_density_kg_per_m_cubed =loaded_point.density

    temperature = loaded_point.temperature
    
    pressure=loaded_point.pressure

    RCC = loaded_point.rcut_couloumb


    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]
    import pickle
    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)

    ### Note:  The electrostatics and Ewald are turned off in the
    # GOMC control file (i.e., False) since the n-alkanes beads in the
    # trappe-ua force field have no charge (i.e., the bead charges are all zero)

    MC_steps=5000
    restart_coor = "${restart_coor}"
    restart_xsc = "${restart_xsc}"
    gomc_control.write_gomc_control_file(charmm, conf_filename='system_npt',  ensemble_type='NPT', RunSteps=MC_steps, Restart=True, \
                                        check_input_files_exist=False, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        Coordinates_box_0="system.pdb",Structure_box_0="system.psf",binCoordinates_box_0=restart_coor,
                                        extendedSystem_box_0=restart_xsc,
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": True,
                                                            "EqSteps" : 1000,
                                                            "AdjSteps":10,
                                                            "Pressure" : float(pressure), 
                                                            "PRNG": int(0),
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": RCC,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.98,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.01,
                                                            "MultiParticleFreq":0.01,
                                                            "OutputName":"system_npt"

                                                                }
                                        )

    """
}

process NPT_equilibration_solvent_system {
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_npt_eq", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint), path(npt_conf), path(pdb), path(psf), path(inp),\
    path("system_nvt_restart.chk"),path("system_nvt_BOX_0_restart.coor"), path("system_nvt_BOX_0_restart.xsc")
    output:
    tuple val(temp_K), val(Rho_kg_per_m_cubed),path(statepoint), path(npt_conf), path(pdb), path(psf), path(inp),\
    path("system_npt_restart.chk"),path("system_npt_BOX_0_restart.coor"), path("system_npt_BOX_0_restart.xsc"), emit: system
    tuple path("NPT_equilibration.log"), path(npt_conf),  emit: record

    shell:
    """
    
    #!/bin/bash
    cat ${npt_conf} > local.conf
    GOMC_CPU_NPT +p${task.cpus} local.conf > NPT_equilibration.log
    """
}

process write_gomc_calibration_confs {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/NVT/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_ewald_calibration_input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed),val(temp_K),path(charmm),val(temp_K),path(statepoint),path(restart_xsc),path(restart_coor)
    output:
    tuple val(Rho_kg_per_m_cubed),val(temp_K),path("ewald_calibration.conf"), emit: ewald_calibration_conf
    script:
    """
    #!/usr/bin/env python
    from typing import List
    from pydantic import BaseModel

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

    loaded_point = load_point_from_json("${statepoint}")

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

    liquid_box_length_Ang = loaded_point.box_length

    liquid_box_density_kg_per_m_cubed =loaded_point.density

    temperature = loaded_point.temperature
    
    pressure=loaded_point.pressure

    RCC = loaded_point.rcut_couloumb

    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]
    import pickle
    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)

    ### Note:  The electrostatics and Ewald are turned off in the
    # GOMC control file (i.e., False) since the n-alkanes beads in the
    # trappe-ua force field have no charge (i.e., the bead charges are all zero)

    MC_steps=10000
    #MC_steps=5000

    restart_coor = "${restart_coor}"
    restart_xsc = "${restart_xsc}"
    if (${params.debugging}):
        MC_steps=2000
    else:
        MC_steps=5000

    #restart_chk = ""
    gomc_control.write_gomc_control_file(charmm, conf_filename='ewald_calibration',  ensemble_type='NVT', RunSteps=MC_steps, Restart=True, \
                                        check_input_files_exist=False, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        Coordinates_box_0="system.pdb",Structure_box_0="system.psf",binCoordinates_box_0=restart_coor,
                                        extendedSystem_box_0=restart_xsc,
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": True,
                                                            "EqSteps" : 1000,
                                                            "AdjSteps":10,
                                                            "Pressure" : float(pressure), 
                                                            "PRNG": int(0),
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": RCC,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.50,
                                                            "RotFreq":0.50,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.00,
                                                            "MultiParticleFreq":0.00,
                                                            "OutputName":"ewald_calibration"

                                                                }
                                        )
    if (${params.debugging}):
        NUM_POINTS = 5.0
    else:
        NUM_POINTS = 50.0
    RCC_START = 10.0
    if (temperature==600):
        percentage = 0.80
    else:
        percentage = 0.80
    RCC_END = (float(liquid_box_length_Ang)/2.0)*percentage
    RCC_DELTA = (RCC_END-RCC_START)/(NUM_POINTS-1)
    ALPHA_START = 0.0
    ALPHA_END = 0.5
    ALPHA_DELTA = (ALPHA_END-ALPHA_START)/(NUM_POINTS-1)
    file1 = open("ewald_calibration.conf", "a")
    #defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file=restart_chk)
    #file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="1000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start=RCC_START,\
    end=RCC_END,delta=RCC_DELTA)
    file1.writelines(defAlphaLine)

    """
}


process GOMC_Ewald_Calibration {
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/NVT/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_ewald_calibration", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path(statepoint),path(pdb), path(psf), path(inp), path(namd_inp),\
    val(temp_K), path(statepoint, stageAs: "copyOfStatepoint.json"), path(restart_xsc), path(restart_coor),\
    val(temp_K), path(ewald_calibration_conf)
    
    output:
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path("Wolf_Calibration_*"), emit: grids
    path("Ewald_Calibration.log"),  emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${ewald_calibration_conf} > local.conf
    GOMC_CPU_NVT +p${task.cpus} local.conf > Ewald_Calibration.log
    """
}

process plot_grids {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/NVT/temperature_${temp_K}_density_${Rho_kg_per_m_cubed}/gomc_ewald_calibration_plots", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), val(temp_K), path(statepoint, stageAs: "copyOfStatepoint.json"), path(restart_xsc), path(restart_coor),val(temp_K), path(grids)
    output:
    tuple path("grids*.png"), path("limited_axes.png"), emit: figs
    tuple val(temp_K), val(Rho_kg_per_m_cubed), path(statepoint), path("convergence_obj.json"), path(restart_xsc), path(restart_coor), emit: convergence
    
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import re
    import numpy as np


    # Function to extract model name from file name using regex
    def extract_model_name(file_name):
        match = re.match(r'Wolf_Calibration_(\\w+)_BOX_\\d+_(\\w+).dat', file_name)
        if match:
            return match.group(1)
        else:
            return None

    # Example list of files
    file_list = "${grids}".split()

    # Determine the number of subplots based on the number of files
    num_subplots = len(file_list)

    # Calculate the number of rows and columns for subplots
    num_rows = (num_subplots + 2) // 3  # Ensure at least 1 row
    num_cols = min(3, num_subplots)

    # Create a new figure with subplots arranged in rows and columns
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows))

    # Create a new figure with subplots arranged in rows and columns
    fig_slopes, axes_slopes = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows))

    fig_convergence, axes_convergence = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows),sharex=True)

    # Flatten the axes array to simplify indexing
    axes = axes.flatten()
    axes_slopes = axes_slopes.flatten()
    axes_convergence = axes_convergence.flatten()
    import random
    extended_markers = [
        'o', 's', '^', 'v', '<', '>', 'D', 'p', '*', 'h', '+', 'x', '|', '_',
        '.', ',', '1', '2', '3', '4', '8', 'H', 'd', 'D', 'P', 'X', 'o', 's', 'p', '*', 'h', '+', 'x'
    ]

    model_dict = {}

    # Iterate over files and plot each DataFrame in a subplot
    for i, file_name in enumerate(file_list):
        random.seed(0)
        # Extract model name
        model_name = extract_model_name(file_name)

        # Load the CSV file into a DataFrame
        df = pd.read_csv(file_name, index_col=0)
        df_slopes = pd.DataFrame(index=df.index, columns=df.columns)

        desired_y_values = np.arange(-2, 2.1, 0.1)

        for col in df.columns:
            x = df.index  # Using DataFrame indices as x-values
            y = df[col]

            # Calculate the slope using numpy's gradient function
            df_slopes[col] = np.gradient(y, x)


        abs_df = df.abs()
        abs_slopes_df = df_slopes.abs()

        normalized_df=(abs_df-abs_df.min())/(abs_df.max()-abs_df.min())
        normalized_slopes_df=(abs_slopes_df-abs_slopes_df.min())/(abs_slopes_df.max()-abs_slopes_df.min())

        # Calculate Euclidean distance for each tuple across all columns
        tuple_df = pd.concat([normalized_df, normalized_slopes_df]).groupby(level=0).apply(lambda x: np.sqrt(np.sum(x**2)))

        # Find the row and column of the minimum entry
        min_entry_location = tuple_df.unstack().idxmin()

        # Extract row and column indices
        min_row, min_col = min_entry_location

        print(f"RCut of Minimum Entry: {min_row}")
        print(f"Alpha of Minimum Entry: {min_col}")
        print(f"Slope of Minimum Entry: ",df_slopes[min_row][min_col])
        print(f"F(alpha) of Minimum Entry: ",df[min_row][min_col])

        # Create a dictionary with the information
        result_dict = {
            'ConvergedRCut': min_row,
            'ConvergedAlpha': min_col,  # Assuming column names are α
            'ConvergedFAlpha': df[min_row][min_col],
            'ConvergedDFDAlpha': df_slopes[min_row][min_col]
        }

        model_dict[model_name]=result_dict

        # Plotting one line per row in the subplot
        #df.plot(ax=axes[i], marker='o', linestyle='-', legend=False)
        # Iterate over columns and plot each with a different marker and color
        for column in df.columns:
            # Generate random color and marker for each column
            color = "#{:06x}".format(random.randint(0, 0xFFFFFF))  # Random hex color
            marker = random.choice(extended_markers)  # Random marker
            fill_style = random.choice(['full', 'left', 'right', 'none'])  # Random fill style

            # Plotting one line per column in the subplot with random color, marker, and fill style
            df[column].plot(
                ax=axes[i],
                marker=marker,
                linestyle='-',
                color=color,
                fillstyle=fill_style,
                legend=False
            )

            # Plotting one line per column in the subplot with random color, marker, and fill style
            df_slopes[column].plot(
                ax=axes_slopes[i],
                marker=marker,
                linestyle='-',
                color=color,
                fillstyle=fill_style,
                legend=False
            )


            # Plotting one line per column in the subplot with random color, marker, and fill style
            axes_convergence[i].plot(
                df[column], 
                df_slopes[column],
                marker=marker,
                linestyle='-',
                color=color,
                fillstyle=fill_style, 
                label=None
            )

        #axes[i].text(result_dict['ConvergedAlpha'], result_dict['ConvergedFAlpha'], "RCut={:.2f}, α={:.2f}".format(float(result_dict['ConvergedRCut']),float(result_dict['ConvergedAlpha'])), fontsize=20, ha='right', va='bottom')
        axes[i].plot(result_dict['ConvergedAlpha'], result_dict['ConvergedFAlpha'], marker='*', markersize=20, linestyle='None', color='red')

        #axes_slopes[i].text(result_dict['ConvergedAlpha'], result_dict['ConvergedDFDAlpha'], "RCut={:.2f}, α={:.2f}".format(float(result_dict['ConvergedRCut']),float(result_dict['ConvergedAlpha'])), fontsize=20, ha='right', va='bottom')
        axes_slopes[i].plot(result_dict['ConvergedAlpha'], result_dict['ConvergedDFDAlpha'], marker='*', markersize=20, linestyle='None', color='red')

        #axes_convergence[i].text(result_dict['ConvergedFAlpha'], result_dict['ConvergedDFDAlpha'], "RCut={:.2f}, α={:.2f}".format(float(result_dict['ConvergedRCut']),float(result_dict['ConvergedAlpha'])), fontsize=20, ha='right', va='bottom')
        axes_convergence[i].plot(result_dict['ConvergedFAlpha'], result_dict['ConvergedDFDAlpha'], marker='*', markersize=20, linestyle='None', color='red')

        # Set plot labels and title 
        axes[i].set_xlabel('α')
        axes[i].set_ylabel('f(α)')
        axes[i].set_title(f'Model: {model_name}')

        # Set plot labels and title
        axes_slopes[i].set_xlabel('α')
        axes_slopes[i].set_ylabel('d/dα f(α)')
        axes_slopes[i].set_title(f'Model: {model_name}')


        # Set plot labels and title
        axes_convergence[i].set_xlabel('f(α)')
        axes_convergence[i].set_ylabel('d/dα f(α)')
        axes_convergence[i].set_title(f'Model: {model_name}')


    # Create a single legend for the entire figure outside the canvas to the right
    handles, labels = axes[0].get_legend_handles_labels()
    #fig.legend(handles, labels, loc='right')
    legend = fig.legend(handles, labels, bbox_to_anchor=(1.10, 0.5), loc='center right')
    # Adjust layout
    fig.tight_layout()
    fig_slopes.tight_layout()
    # Save the figure as 'grids.png'
    fig.savefig('grids.png', bbox_inches='tight', bbox_extra_artists=[legend])
    fig_slopes.savefig('grids_slopes.png', bbox_inches='tight', bbox_extra_artists=[legend])
    fig_convergence.savefig('grids_convergence.png', bbox_inches='tight', bbox_extra_artists=[legend])

    # Show the plot

    # Create a second figure with subplots and limited y-axes to -1 to +1
    fig2, axes2 = fig, axes

    # Flatten the axes array to simplify indexing
    axes2 = axes2.flatten()

    # Iterate over subplots in the second figure and set y-axis limits
    for ax in axes2:
        ax.set_ylim(-1, 1)

    # Save the second figure as 'limited_axes.png'
    plt.savefig('limited_axes.png', bbox_inches='tight')



    # Show the plots
    plt.show()

    from typing import Dict, Union
    from pydantic import BaseModel, Field
    import json

    class InnerModel(BaseModel):
        ConvergedRCut: float
        ConvergedAlpha: float
        ConvergedFAlpha: float
        ConvergedDFDAlpha: float

    class FooBarModel(BaseModel):
        models: Dict[str, InnerModel]

    # Create an instance of the model
    convergence_obj = FooBarModel(
        models=model_dict
    )

    # Serialize the Pydantic object to JSON
    with open("convergence_obj.json", 'w') as file:
        file.write(convergence_obj.model_dump_json())
    """
}


process plot_grids_two_box {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/calibration/plots", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), path(grids)
    output:
    tuple path("grids*.png"), path("limited_axes.png"), emit: figs
    tuple val(temp_K), path("convergence_obj.json"), emit: convergence
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import re
    import numpy as np
    #rel_error_weight = 0.75
    #slope_weight = 0.25

    rel_error_weight = 1.0
    slope_weight = 0.0

    # Function to extract model name from file name using regex
    def extract_model_name(file_name):
        match = re.match(r'Wolf_Calibration_(\\w+)_BOX_(\\d+)_(.*)\\.dat', file_name)
        if match:
            return match.group(1)
        else:
            return None

    # Function to extract model name from file name using regex
    def extract_box(file_name):
        match = re.match(r'Wolf_Calibration_(\\w+)_BOX_(\\d+)_(.*)\\.dat', file_name)
        if match:
            return match.group(2)
        else:
            return None

    # Example list of files
    file_list = "${grids}".split()

    # Determine the number of subplots based on the number of files
    num_subplots = len(file_list)

    # Calculate the number of rows and columns for subplots
    num_rows = (num_subplots + 1) // 2  # Ensure at least 1 row, with 2 columns per row
    num_cols = min(2, num_subplots)  # 2 columns per row

    # Create a new figure with subplots arranged in rows and columns
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows))

    # Create a new figure with subplots arranged in rows and columns
    fig_slopes, axes_slopes = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows))

    fig_convergence, axes_convergence = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows), sharex=True)

    # Flatten the axes array to simplify indexing
    axes = axes.flatten()
    axes_slopes = axes_slopes.flatten()
    axes_convergence = axes_convergence.flatten()

    # Hide any extra subplots
    for i in range(num_subplots, num_rows * num_cols):
        fig.delaxes(axes[i])
        fig_slopes.delaxes(axes_slopes[i])
        fig_convergence.delaxes(axes_convergence[i])

    import random
    extended_markers = [
        'o', 's', '^', 'v', '<', '>', 'D', 'p', '*', 'h', '+', 'x', '|', '_',
        '.', ',', '1', '2', '3', '4', '8', 'H', 'd', 'D', 'P', 'X', 'o', 's', 'p', '*', 'h', '+', 'x'
    ]

    model_dict = {}

    # Iterate over files and plot each DataFrame in a subplot
    for i, file_name in enumerate(file_list):
        random.seed(0)
        # Extract model name
        model_name = extract_model_name(file_name)
        box = extract_box(file_name)

        # Load the CSV file into a DataFrame
        df = pd.read_csv(file_name, index_col=0)
        df_slopes = pd.DataFrame(index=df.index, columns=df.columns)

        desired_y_values = np.arange(-2, 2.1, 0.1)

        for col in df.columns:
            x = df.index  # Using DataFrame indices as x-values
            y = df[col]

            # Calculate the slope using numpy's gradient function
            df_slopes[col] = np.gradient(y, x)


        abs_df = df.abs()
        abs_slopes_df = df_slopes.abs()

        normalized_df=(abs_df-abs_df.min())/(abs_df.max()-abs_df.min())
        normalized_slopes_df=(abs_slopes_df-abs_slopes_df.min())/(abs_slopes_df.max()-abs_slopes_df.min())

        normalized_df = normalized_df.multiply(rel_error_weight)
        normalized_slopes_df = normalized_slopes_df.multiply(slope_weight)

        # Calculate Euclidean distance for each tuple across all columns
        tuple_df = pd.concat([normalized_df, normalized_slopes_df]).groupby(level=0).apply(lambda x: np.sqrt(np.sum(x**2)))

        # Find the row and column of the minimum entry
        min_entry_location = tuple_df.unstack().idxmin()

        # Extract row and column indices
        min_row, min_col = min_entry_location
        print(f"Model name: {model_name}")
        print(f"Box: {box}")
        print(f"RCut of Minimum Entry: {min_row}")
        print(f"Alpha of Minimum Entry: {min_col}")
        print(f"Slope of Minimum Entry: ",df_slopes[min_row][min_col])
        print(f"F(alpha) of Minimum Entry: ",df[min_row][min_col])

        # Create a dictionary with the information
        result_dict = {
            'ConvergedRCut': min_row,
            'ConvergedAlpha': min_col,  # Assuming column names are α
            'ConvergedFAlpha': df[min_row][min_col],
            'ConvergedDFDAlpha': df_slopes[min_row][min_col]
        }
        # Add the result_dict to the model_dict under the appropriate model_name and each box_value
        if model_name not in model_dict:
            model_dict[model_name] = {}

        #for box_value in box_values:
        #    model_dict[model_name][box_value] = result_dict
        model_dict[model_name][box]=result_dict

        # Plotting one line per row in the subplot
        #df.plot(ax=axes[i], marker='o', linestyle='-', legend=False)
        # Iterate over columns and plot each with a different marker and color
        for column in df.columns:
            # Generate random color and marker for each column
            color = "#{:06x}".format(random.randint(0, 0xFFFFFF))  # Random hex color
            marker = random.choice(extended_markers)  # Random marker
            fill_style = random.choice(['full', 'left', 'right', 'none'])  # Random fill style

            # Plotting one line per column in the subplot with random color, marker, and fill style
            df[column].plot(
                ax=axes[i],
                marker=marker,
                linestyle='-',
                color=color,
                fillstyle=fill_style,
                legend=False
            )

            # Plotting one line per column in the subplot with random color, marker, and fill style
            df_slopes[column].plot(
                ax=axes_slopes[i],
                marker=marker,
                linestyle='-',
                color=color,
                fillstyle=fill_style,
                legend=False
            )


            # Plotting one line per column in the subplot with random color, marker, and fill style
            axes_convergence[i].plot(
                df[column], 
                df_slopes[column],
                marker=marker,
                linestyle='-',
                color=color,
                fillstyle=fill_style, 
                label=None
            )

        #axes_slopes[i].text(result_dict['ConvergedAlpha'], result_dict['ConvergedDFDAlpha'], "RCut={:.2f}, α={:.2f}".format(float(result_dict['ConvergedRCut']),float(result_dict['ConvergedAlpha'])), fontsize=20, ha='right', va='bottom')
        axes_slopes[i].plot(result_dict['ConvergedAlpha'], result_dict['ConvergedDFDAlpha'], marker='*', markersize=20, linestyle='None', color='red')

        #axes_convergence[i].text(result_dict['ConvergedFAlpha'], result_dict['ConvergedDFDAlpha'], "RCut={:.2f}, α={:.2f}".format(float(result_dict['ConvergedRCut']),float(result_dict['ConvergedAlpha'])), fontsize=20, ha='right', va='bottom')
        axes_convergence[i].plot(result_dict['ConvergedFAlpha'], result_dict['ConvergedDFDAlpha'], marker='*', markersize=20, linestyle='None', color='red')

        #axes[i].text(result_dict['ConvergedAlpha'], result_dict['ConvergedFAlpha'], "RCut={:.2f}, α={:.2f}".format(float(result_dict['ConvergedRCut']),float(result_dict['ConvergedAlpha'])), fontsize=20, ha='right', va='bottom')
        axes[i].plot(result_dict['ConvergedAlpha'], result_dict['ConvergedFAlpha'], marker='*', markersize=20, linestyle='None', color='red')

        # Set plot labels and title
        axes_slopes[i].set_xlabel('α')
        axes_slopes[i].set_ylabel('d/dα f(α)')
        axes_slopes[i].set_title(f'Model: {model_name}; Box: {box}')

        # Set plot labels and title
        axes_convergence[i].set_xlabel('f(α)')
        axes_convergence[i].set_ylabel('d/dα f(α)')
        axes_convergence[i].set_title(f'Model: {model_name}; Box: {box}')

        # Set plot labels and title 
        axes[i].set_xlabel('α')
        axes[i].set_ylabel('f(α)')
        axes[i].set_title(f'Model: {model_name}; Box: {box}')


    # Create a single legend for the entire figure outside the canvas to the right
    handles, labels = axes[0].get_legend_handles_labels()
    #fig.legend(handles, labels, loc='right')
    legend = fig.legend(handles, labels, bbox_to_anchor=(1.10, 0.5), loc='center right')
    # Adjust layout
    fig.tight_layout()
    fig_slopes.tight_layout()
    # Save the figure as 'grids.png'
    fig.savefig('grids.png', bbox_inches='tight', bbox_extra_artists=[legend])
    fig_slopes.savefig('grids_slopes.png', bbox_inches='tight', bbox_extra_artists=[legend])
    fig_convergence.savefig('grids_convergence.png', bbox_inches='tight', bbox_extra_artists=[legend])

    # Show the plot

    # Create a second figure with subplots and limited y-axes to -1 to +1
    fig2, axes2 = fig, axes

    # Flatten the axes array to simplify indexing
    axes2 = axes2.flatten()

    # Iterate over subplots in the second figure and set y-axis limits
    for ax in axes2:
        ax.set_ylim(-2, 2)

    # Save the second figure as 'limited_axes.png'
    fig2.savefig('limited_axes.png', bbox_inches='tight')

    # Show the plots
    plt.show()

    from typing import Dict, Union
    from pydantic import BaseModel, Field
    import json

    class InnerModel(BaseModel):
        ConvergedRCut: float
        ConvergedAlpha: float
        ConvergedFAlpha: float
        ConvergedDFDAlpha: float

    class FooBarModel(BaseModel):
        models: Dict[str, Dict[int, InnerModel]]  # Adjusted for the nested dictionary

    # Convert the result_dict to InnerModel object
    inner_result_dict = InnerModel(**result_dict)

    # Create an instance of the model
    convergence_obj = FooBarModel(
        models=model_dict
    )

    # Serialize the Pydantic object to JSON
    with open("convergence_obj.json", 'w') as file:
        file.write(convergence_obj.json())
    """
}



process write_gemc_ewald_confs {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/ewald/input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2),\
    path(charmm),path(statepoint1, stageAs: "statepoint1.json"), path(statepoint2, stageAs: "statepoint2.json")    
    output:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2),\
    path("in_GEMC_NVT.conf"), emit:gemc_conf
    tuple val(temp_K), path(charmm),path(statepoint1), path(statepoint2), emit:charmm
    script:
    """
    #!/usr/bin/env python
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

    import pickle
    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)


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
        gomc_steps_equilibration = 10000000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 10000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = 100000 # Monte Carlo Steps between console output
    coordinate_output_freq = 100000 # # set value for paper = 50 * 10**3
    restart_output_freq = 100000 # # set value for paper = 50 * 10**3
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
    RcutCoulomb_box_0=loaded_point1.rcut_couloumb
    RcutCoulomb_box_1=loaded_point2.rcut_couloumb

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
                                                            "RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            "RcutCoulomb_box_1": RcutCoulomb_box_1,
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
                                                            "RestartFreq": [True, restart_output_freq],
                                                            "CheckpointFreq": [True, restart_output_freq],
                                                            "DCDFreq": [True, coordinate_output_freq],
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


process write_gemc_ewald_confs_calibrate {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/calibrate/confs", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2),\
    path(charmm),path(statepoint1, stageAs: "statepoint1.json"), path(statepoint2, stageAs: "statepoint2.json")    
    output:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2),\
    path("in_GEMC_NVT.conf"), emit:gemc_conf
    tuple val(temp_K), path(charmm),path(statepoint1), path(statepoint2), emit:charmm
    script:
    """
    #!/usr/bin/env python
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

    import pickle
    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)


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
        gomc_steps_equilibration = 10000000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 10000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = 100000 # Monte Carlo Steps between console output
    coordinate_output_freq = 100000 # # set value for paper = 50 * 10**3
    restart_output_freq = 100000 # # set value for paper = 50 * 10**3
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
    RcutCoulomb_box_0=loaded_point1.rcut_couloumb
    RcutCoulomb_box_1=loaded_point2.rcut_couloumb

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
                                                            "RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            "RcutCoulomb_box_1": RcutCoulomb_box_1,
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
                                                            "RestartFreq": [True, restart_output_freq],
                                                            "CheckpointFreq": [True, restart_output_freq],
                                                            "DCDFreq": [True, coordinate_output_freq],
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
    if (${params.debugging}):
        NUM_POINTS = 5.0
    else:
        NUM_POINTS = 50.0
    percentage = 0.80
    ALPHA_START = 0.0
    ALPHA_END = 0.5
    ALPHA_DELTA = (ALPHA_END-ALPHA_START)/(NUM_POINTS-1)
    RCC_START = 10.0
    RCC_END_BOX_0 = (float(liquid_box_0_length_Ang)/2.0)*percentage
    RCC_DELTA_BOX_0 = (RCC_END_BOX_0-RCC_START)/(NUM_POINTS-1)
    RCC_END_BOX_1 = (float(liquid_box_1_length_Ang)/2.0)*percentage
    RCC_DELTA_BOX_1 = (RCC_END_BOX_1-RCC_START)/(NUM_POINTS-1)
    file1 = open("calibration.conf", "a")
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="1000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start=RCC_START,\
    end=RCC_END_BOX_0,delta=RCC_DELTA_BOX_0)
    file1.writelines(defAlphaLine)

    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="1",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="1",start=RCC_START,\
    end=RCC_END_BOX_1,delta=RCC_DELTA_BOX_1)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file="${chk}")
    #file1.writelines(defAlphaLine)

    charmm.write_inp()
    """
}


process write_gemc_ewald_calibration_confs {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/methods/calibration/input", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input:
    tuple val(temp_K),\
    path(rstpdb1), path(rstpsf1), path(rstpdb2), path(rstpsf2),\
    path(rstcoor1),path(rstxsc1), path(rstcoor2), path(rstxsc2), path(chk),\
    path(charmm),path(statepoint1), path(statepoint2)
    output:
    tuple val(temp_K),\
    path(rstpdb1), path(rstpsf1), path(rstpdb2), path(rstpsf2),\
    path(rstcoor1),path(rstxsc1), path(rstcoor2), path(rstxsc2), path(chk),\
    path("system.inp"), path("calibration.conf"), emit:gemc_calibration

    script:
    """
    #!/usr/bin/env python
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

    liquid_box_0_length_Ang = loaded_point1.box_length
    liquid_box_1_length_Ang = loaded_point2.box_length

    import pickle
    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)


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
        gomc_steps_equilibration = 10000000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 10000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = 100000 # Monte Carlo Steps between console output
    coordinate_output_freq = 100000 # # set value for paper = 50 * 10**3
    restart_output_freq = 100000 # # set value for paper = 50 * 10**3
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
    RcutCoulomb_box_0=loaded_point1.rcut_couloumb
    RcutCoulomb_box_1=loaded_point2.rcut_couloumb

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
    output_file_prefix="GOMC_GEMC_Calibration"


    gomc_control.write_gomc_control_file(charmm, 'calibration.conf', 'GEMC_NVT', MC_steps, ${temp_K},
                                        Restart=True,
                                        check_input_files_exist=False,
                                        Coordinates_box_0="${rstpdb1}",
                                        Structure_box_0="${rstpsf1}",
                                        binCoordinates_box_0="${rstcoor1}",
                                        extendedSystem_box_0="${rstxsc1}",
                                        Coordinates_box_1="${rstpdb2}",
                                        Structure_box_1="${rstpsf2}",
                                        binCoordinates_box_1="${rstcoor2}",
                                        extendedSystem_box_1="${rstxsc2}",
                                        input_variables_dict={"VDWGeometricSigma": False,
                                                            "Ewald": True,
                                                            "ElectroStatic": True,
                                                            "PRNG": int(0),
                                                            "Pressure": None,
                                                            "Potential": "VDW",
                                                            "LRC": LRC,
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            "RcutCoulomb_box_1": RcutCoulomb_box_1,
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
                                                            "RestartFreq": [True, restart_output_freq],
                                                            "CheckpointFreq": [True, restart_output_freq],
                                                            "DCDFreq": [True, coordinate_output_freq],
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

    if (${params.debugging}):
        NUM_POINTS = 5.0
    else:
        NUM_POINTS = 50.0
    percentage = 0.80
    ALPHA_START = 0.0
    ALPHA_END = 0.5
    ALPHA_DELTA = (ALPHA_END-ALPHA_START)/(NUM_POINTS-1)
    RCC_START = 10.0
    RCC_END_BOX_0 = (float(liquid_box_0_length_Ang)/2.0)*percentage
    RCC_DELTA_BOX_0 = (RCC_END_BOX_0-RCC_START)/(NUM_POINTS-1)
    RCC_END_BOX_1 = (float(liquid_box_1_length_Ang)/2.0)*percentage
    RCC_DELTA_BOX_1 = (RCC_END_BOX_1-RCC_START)/(NUM_POINTS-1)
    file1 = open("calibration.conf", "a")
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="1000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start=RCC_START,\
    end=RCC_END_BOX_0,delta=RCC_DELTA_BOX_0)
    file1.writelines(defAlphaLine)

    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="1",start=ALPHA_START,\
    end=ALPHA_END,delta=ALPHA_DELTA)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="1",start=RCC_START,\
    end=RCC_END_BOX_1,delta=RCC_DELTA_BOX_1)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file="${chk}")
    #file1.writelines(defAlphaLine)

    charmm.write_inp()
    """
}


process write_gemc_production_confs {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/methods/${METHOD}/input", mode: 'copy', overwrite: false
    cpus 1

    debug false
    input:
    tuple val(temp_K), path(charmm),path(xsc1), path(coor1), path(xsc2), path(coor2),\
    path(convergenceJSON1, stageAs: "convergenceJSON1.json"), path(convergenceJSON2, stageAs: "convergenceJSON2.json"),\
    val(METHOD)

    output:
    tuple val(temp_K), val(METHOD),path("in_GEMC_NVT.conf"),emit:gemc_conf

    script:
    """
    #!/usr/bin/env python
    from typing import List
    from pydantic import BaseModel
    print("hello from ", $temp_K, "$METHOD" )

    from typing import Dict, Union
    from pydantic import BaseModel, Field
    import json

    class InnerModel(BaseModel):
        ConvergedRCut: float
        ConvergedAlpha: float
        ConvergedFAlpha: float
        ConvergedDFDAlpha: float

    class FooBarModel(BaseModel):
        models: Dict[str, InnerModel]


    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> FooBarModel:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return FooBarModel.model_validate_json(json_data)
    
    loaded_point1 = load_point_from_json("${convergenceJSON1}")
    loaded_point2 = load_point_from_json("${convergenceJSON2}")
    print(loaded_point1)
    print(loaded_point2)

    import pickle
    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)


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
        gomc_steps_equilibration = 10000000
    gomc_steps_production = gomc_steps_equilibration # set value for paper = 1 * 10**6
    console_output_freq = 100 # Monte Carlo Steps between console output
    pressure_calc_freq = 1000 # Monte Carlo Steps for pressure calculation
    block_ave_output_freq = 100000 # Monte Carlo Steps between console output
    coordinate_output_freq = 100 # # set value for paper = 50 * 10**3
    restart_output_freq = 100000 # # set value for paper = 50 * 10**3
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
    RcutCoulomb_box_0=loaded_point1.models["${METHOD}"].ConvergedRCut
    RcutCoulomb_box_1=loaded_point2.models["${METHOD}"].ConvergedRCut

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
                                                            "Ewald": False,
                                                            "ElectroStatic": True,
                                                            "PRNG": int(0),
                                                            "Pressure": None,
                                                            "Potential": "VDW",
                                                            "LRC": LRC,
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": RcutCoulomb_box_0,
                                                            "RcutCoulomb_box_1": RcutCoulomb_box_1,
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
                                                            "RestartFreq": [True, restart_output_freq],
                                                            "CheckpointFreq": [True, restart_output_freq],
                                                            "DCDFreq": [True, coordinate_output_freq],
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

    kind_pot = "${METHOD}".split('_')

    file1 = open("in_GEMC_NVT.conf", "a")
    defAlphaLine = "{box}\\t{val}\\n".format(box="Wolf", val="True")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\n".format(box="WolfKind", val=kind_pot[0])
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\n".format(box="WolfPotential", val=kind_pot[1])
    file1.writelines(defAlphaLine)
    defAlphaLine = "{key}\\t{box}\\t{val}\\n".format(key="WolfAlpha",box="0", val=loaded_point1.models["${METHOD}"].ConvergedAlpha)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{key}\\t{box}\\t{val}\\n".format(key="WolfAlpha",box="1", val=loaded_point2.models["${METHOD}"].ConvergedAlpha)
    file1.writelines(defAlphaLine)

    """
}

process GOMC_GEMC_Production {
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/wolf/methods/${METHOD}/production", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K), val(METHOD), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp),\
    path(xsc1),path(coor1),path(xsc2),path(coor2),path(chk),path(gemc_conf)
    output:
    path("*"), emit: grids
    tuple val(temp_K),val(METHOD),path("GOMC_GEMC_Production.log"),  emit: record

    shell:
    """
    
    #!/bin/bash
    cat ${gemc_conf} > local.conf
    GOMC_CPU_GEMC +p${task.cpus} local.conf > GOMC_GEMC_Production.log
    """
}


process GOMC_GEMC_Production_Replica {
    //cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/ewald/production", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K), path(pdb1), path(psf1), path(pdb2), path(psf2), path(inp), \
    path(xsc1), path(coor1), path(xsc2), path(coor2), path(gemc_conf)
    output:
    //path("*"), emit: grids
    tuple val(temp_K),\
    path("GOMC_GEMC_Production_BOX_0_restart.xsc"),\
    path("GOMC_GEMC_Production_BOX_1_restart.xsc"),\
    path("GOMC_GEMC_Production_BOX_0_restart.coor"),\
    path("GOMC_GEMC_Production_BOX_1_restart.coor"),\
    path("GOMC_GEMC_Production_BOX_0_restart.pdb"),\
    path("GOMC_GEMC_Production_BOX_1_restart.pdb"),\
    path("GOMC_GEMC_Production_BOX_0_restart.psf"),\
    path("GOMC_GEMC_Production_BOX_1_restart.psf"),\
    path("GOMC_GEMC_Production_restart.chk"),  emit: restart_files
    tuple val(temp_K), val("EWALD"), path("GOMC_GEMC_Production.log"),  emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${gemc_conf} > local.conf
    GOMC_CPU_GEMC +p${task.cpus} local.conf > GOMC_GEMC_Production.log
    """
}


process GOMC_GEMC_Calibration {
    cache 'lenient'
    fair true
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/calibration/production", mode: 'copy', overwrite: false
    cpus 8
    debug false
    input:
    tuple val(temp_K),path(pdb1), path(psf1), path(pdb2), path(psf2),path(inp),\
    path(xsc1),path(coor1),path(xsc2),path(coor2),path(chk),path(gemc_conf)
    output:
    tuple val(temp_K), path("Wolf_Calibration_*"), emit: grids
    tuple val(temp_K),path("GOMC_GEMC_Production.log"),  emit: record
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
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/wolf/methods/${METHOD}/extracted/density", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input:
    tuple val(temp_K),val(METHOD),path(logfile)
    output:
    tuple path("${temp_K}_${METHOD}_BOX_0.csv"),path("${temp_K}_${METHOD}_BOX_1.csv"),\
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
    df = pd.DataFrame(densities_box_0, index=None, columns=['${temp_K}_${METHOD}_BOX_0'])
    df.to_csv("${temp_K}_${METHOD}_BOX_0.csv", header=True, sep=' ', index=False)

    steps_box_1, densities_box_1 = extract("$logfile","STAT_1")
    df = pd.DataFrame(densities_box_1, index=None, columns=['${temp_K}_${METHOD}_BOX_1'])
    df.to_csv("${temp_K}_${METHOD}_BOX_1.csv", header=True, sep=' ', index=False)
    """
}


process Collate_GOMC_GEMC_Production { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/density", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path("*")
    output: path("merged_data.csv")

    shell:
    """
    paste * > merged_data.csv
    """
}


process Extract_Vapor_Pressure_GOMC_GEMC_Production {
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/GEMC/temperature_${temp_K}_gemc/wolf/methods/${METHOD}/extracted/vapor_pressure", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input:
    tuple val(temp_K),val(METHOD),path(logfile)
    output:
    tuple path("${temp_K}_${METHOD}_BOX_0.csv"),path("${temp_K}_${METHOD}_BOX_1.csv"),\
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
                            densities.append(float(line.split()[3]))
                        except:
                            print(line)
                            print("An exception occurred") 
        except:
            print("Cant open",filename) 
        steps_np = np.array(steps)
        densities_np = np.array(densities)
        return steps_np, densities_np

    steps_box_0, densities_box_0 = extract("$logfile","STAT_0")
    df = pd.DataFrame(densities_box_0, index=None, columns=['${temp_K}_${METHOD}_BOX_0'])
    df.to_csv("${temp_K}_${METHOD}_BOX_0.csv", header=True, sep=' ', index=False)

    steps_box_1, densities_box_1 = extract("$logfile","STAT_1")
    df = pd.DataFrame(densities_box_1, index=None, columns=['${temp_K}_${METHOD}_BOX_1'])
    df.to_csv("${temp_K}_${METHOD}_BOX_1.csv", header=True, sep=' ', index=False)
    """
}


process Collate_GOMC_GEMC_Production_VP { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/vapor_pressure", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path("*")
    output: path("merged_data.csv")

    shell:
    """
    paste * > merged_data.csv
    """
}



process Plot_GOMC_GEMC_Production { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/density", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: tuple path ("box_0.png"),path("box_1.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h'],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none'],  # Alternating fill pattern
        ['red', 'black', 'blue', 'green', 'purple', 'purple', 'orange', 'green', 'red', 'red'],  # Alternating colors
    )}

    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        if 'BOX_0' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'blue'))  # Default values
            plt.plot(df.index, df[col], label=f'{col}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

            # Plot a solid horizontal line at the initial point of each temperature
            initial_value = df.loc[0, col]
            plt.axhline(y=initial_value, color="black", linestyle='--', linewidth=2)


    plt.xlabel('MC Step', fontsize=18)
    plt.ylabel('Density', fontsize=18)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Box 0 Density/MC Step', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('box_0.png', bbox_inches='tight')

    # Display the plot
    plt.show()

    # Plot for Box 1
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 1
    for idx, col in enumerate(df.columns):
        if 'BOX_1' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'blue'))  # Default values
            plt.plot(df.index, df[col], label=f'{col}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

            # Plot a solid horizontal line at the initial point of each temperature
            initial_value = df.loc[0, col]
            plt.axhline(y=initial_value, color="black", linestyle='--', linewidth=2)


    plt.xlabel('MC Step', fontsize=18)
    plt.ylabel('Density', fontsize=18)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Box 1 Density/MC Step', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('box_1.png', bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}


process Plot_GOMC_GEMC_Production_Ewald { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/ewald/density", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: tuple path ("box_0.png"),path("box_1.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h'],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none'],  # Alternating fill pattern
        ['red', 'black', 'blue', 'green', 'purple', 'purple', 'orange', 'green', 'red', 'red'],  # Alternating colors
    )}

    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        if 'BOX_0' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'blue'))  # Default values
            plt.plot(df.index, df[col], label=f'{col}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

            # Plot a solid horizontal line at the initial point of each temperature
            initial_value = df.loc[0, col]
            plt.axhline(y=initial_value, color="black", linestyle='--', linewidth=2)


    plt.xlabel('MC Step', fontsize=18)
    plt.ylabel('Density', fontsize=18)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Box 0 Density/MC Step', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('box_0.png', bbox_inches='tight')

    # Display the plot
    plt.show()

    # Plot for Box 1
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 1
    for idx, col in enumerate(df.columns):
        if 'BOX_1' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'blue'))  # Default values
            plt.plot(df.index, df[col], label=f'{col}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

            # Plot a solid horizontal line at the initial point of each temperature
            initial_value = df.loc[0, col]
            plt.axhline(y=initial_value, color="black", linestyle='--', linewidth=2)


    plt.xlabel('MC Step', fontsize=18)
    plt.ylabel('Density', fontsize=18)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Box 1 Density/MC Step', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('box_1.png', bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}


process Plot_GOMC_GEMC_Production_VLE { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/density", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: path ("vle.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h','s',],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':','-'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full'],  # Alternating fill pattern
        ['red', 'red', 'blue', 'blue', 'green', 'yellow', 'orange', 'green', 'green', 'grey', 'purple'],  # Alternating colors
    )}
    method_data = {}
    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        temperature = col.split('_')[0]
        method = col.split('_')[1] + '_' + col.split('_')[2]
        # Append temperature and density to the corresponding method in the dictionary
        if method not in method_data:
            method_data[method] = {'temperature': [], 'density': []}
        method_data[method]['temperature'].append(int(temperature))
        method_data[method]['density'].append(df[col].mean())
        if (method == "RAHBARI_DSF"):
            method = "EWALD"
            if method not in method_data:
                method_data[method] = {'temperature': [], 'density': []}
            method_data[method]['temperature'].append(int(temperature))
            method_data[method]['density'].append(df[col][0])

    # Sort temperature and density lists for each method by density
    for method, data in method_data.items():
        temp_density = list(zip(data['temperature'], data['density']))
        temp_density.sort(key=lambda x: x[1])  # Sort by density
        method_data[method]['temperature'], method_data[method]['density'] = zip(*temp_density)

    # Plotting sorted data
    plt.plot(method_data["EWALD"]['density'], method_data["EWALD"]['temperature'], label=f'EWALD', linestyle='-', color="black", marker='o', fillstyle="full")
    for method, data in method_data.items():
        if method != "EWALD":
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'black'))
            plt.plot(data['density'], data['temperature'], label=f'{method}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)


    plt.xlabel('Density', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    lgd = plt.legend(by_label.values(), by_label.keys(),fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('VLE', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('vle.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}


process Plot_GOMC_GEMC_Production_VLE_Per_Density { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/density", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: path ("per_density.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h','s',],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':','-'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full'],  # Alternating fill pattern
        ['red', 'red', 'blue', 'blue', 'green', 'yellow', 'orange', 'green', 'green', 'grey', 'purple'],  # Alternating colors
    )}
    method_data = {}
    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        temperature = col.split('_')[0]
        method = col.split('_')[1] + '_' + col.split('_')[2]
        # Append temperature and density to the corresponding method in the dictionary
        if method not in method_data:
            method_data[method] = {'temperature': [], 'density': []}
        method_data[method]['temperature'].append(int(temperature))
        method_data[method]['density'].append(df[col].mean())
        if (method == "RAHBARI_DSF"):
            method = "EWALD"
            if method not in method_data:
                method_data[method] = {'temperature': [], 'density': []}
            method_data[method]['temperature'].append(int(temperature))
            method_data[method]['density'].append(df[col][0])

    # Sort temperature and density lists for each method by density
    for method, data in method_data.items():
        temp_density = list(zip(data['temperature'], data['density']))
        temp_density.sort(key=lambda x: x[1])  # Sort by density
        method_data[method]['temperature'], method_data[method]['density'] = zip(*temp_density)

    import matplotlib.pyplot as plt

    # Get unique density values
    unique_densities = sorted(set(method_data["EWALD"]["density"]))

    # Create subplots based on the number of unique density values
    num_subplots = len(unique_densities)
    fig, axs = plt.subplots(1, num_subplots, figsize=(15, 5))

    markers = ['o', 's', '^', 'D', 'v', '<', ]
    colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple'],  # Alternating colors

    # Iterate over each subplot
    # Iterate over each subplot
    for i, density in enumerate(unique_densities):
        axs[i].set_title(f'Ewald Density: {density}')
        axs[i].set_ylabel('Density')
        axs[i].tick_params(axis='x', rotation=90)  # Adjust rotation angle as needed
        # Get temperature and method values for the current density
        methods = []
        densities = []
        ew_dens = 0
        for method, data in method_data.items():
            if method != "EWALD":
                methods.append(method)
                densities.append(data['density'][i])
            else:
                ew_dens = data['density'][i]
        # Plotting bar chart for each method
        #axs[i].bar(methods, densities, color='blue')
        axs[i].scatter(methods, densities)
        axs[i].axhline(y=ew_dens, color="black", linestyle='--', linewidth=2)
    # Save the plot as a PNG file
    plt.savefig('per_density.png', bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}


process Plot_GOMC_GEMC_Production_VP { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/vapor_pressure", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: path("vapor_pressure.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h'],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none'],  # Alternating fill pattern
        ['red', 'black', 'blue', 'green', 'purple', 'purple', 'orange', 'green', 'red', 'red'],  # Alternating colors
    )}

    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        if 'BOX_0' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'blue'))  # Default values
            plt.plot(df.index, df[col], label=f'{col}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

            # Plot a solid horizontal line at the initial point of each temperature
            initial_value = df.loc[0, col]
            plt.axhline(y=initial_value, color="black", linestyle='--', linewidth=2)


    plt.xlabel('MC Step', fontsize=18)
    plt.ylabel('Vapor Pressure', fontsize=18)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Box 0 Density/MC Step', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('box_0.png', bbox_inches='tight')

    # Display the plot
    plt.show()

    # Plot for Box 1
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 1
    for idx, col in enumerate(df.columns):
        if 'BOX_1' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'blue'))  # Default values
            plt.plot(df.index, df[col], label=f'{col}', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

            # Plot a solid horizontal line at the initial point of each temperature
            initial_value = df.loc[0, col]
            plt.axhline(y=initial_value, color="black", linestyle='--', linewidth=2)


    plt.xlabel('MC Step', fontsize=18)
    plt.ylabel('Density', fontsize=18)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Box 1 Density/MC Step', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('box_1.png', bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}


process Plot_GOMC_GEMC_Production_VLE_VP { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/vapor_pressure", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: path ("vapor_pressure.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h','s',],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':','-'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full'],  # Alternating fill pattern
        ['red', 'red', 'blue', 'blue', 'green', 'yellow', 'orange', 'green', 'green', 'grey', 'purple'],  # Alternating colors
    )}
    method_data_box_0 = {}
    method_data_box_1 = {}
    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        if 'BOX_0' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            # Append temperature and density to the corresponding method in the dictionary
            if method not in method_data_box_0:
                method_data_box_0[method] = {'temperature': [], 'vapor_pressure': []}
            method_data_box_0[method]['temperature'].append(int(temperature))
            method_data_box_0[method]['vapor_pressure'].append(df[col].mean())
            if (method == "RAHBARI_DSF"):
                method = "EWALD"
                if method not in method_data_box_0:
                    method_data_box_0[method] = {'temperature': [], 'vapor_pressure': []}
                method_data_box_0[method]['temperature'].append(int(temperature))
                method_data_box_0[method]['vapor_pressure'].append(df[col][0])


        if 'BOX_1' in col:
            temperature = col.split('_')[0]
            method = col.split('_')[1] + '_' + col.split('_')[2]
            # Append temperature and density to the corresponding method in the dictionary
            if method not in method_data_box_1:
                method_data_box_1[method] = {'temperature': [], 'vapor_pressure': []}
            method_data_box_1[method]['temperature'].append(int(temperature))
            method_data_box_1[method]['vapor_pressure'].append(df[col].mean())
            if (method == "RAHBARI_DSF"):
                method = "EWALD"
                if method not in method_data_box_1:
                    method_data_box_1[method] = {'temperature': [], 'vapor_pressure': []}
                method_data_box_1[method]['temperature'].append(int(temperature))
                method_data_box_1[method]['vapor_pressure'].append(df[col][0])

    # Plotting sorted data
    plt.plot(method_data_box_0["EWALD"]['temperature'], method_data_box_0["EWALD"]['vapor_pressure'], label=f'EWALD Liquid', linestyle='-', color="black", marker='o', fillstyle="full")
    for method, data in method_data_box_0.items():
        if method != "EWALD":
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'black'))
            plt.plot(data['temperature'], data['vapor_pressure'],  label=f'{method} Liquid', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

    plt.plot(method_data_box_1["EWALD"]['temperature'], method_data_box_1["EWALD"]['vapor_pressure'], label=f'EWALD Vapor', linestyle='-', color="black", marker='o', fillstyle="full")
    for method, data in method_data_box_1.items():
        if method != "EWALD":
            marker, linestyle, fillstyle, color = method_properties.get(method, ('o', '-', 'full', 'black'))
            plt.plot(data['temperature'], data['vapor_pressure'],  label=f'{method} Vapor', linestyle=linestyle, color=color, marker=marker, fillstyle=fillstyle)

    plt.xlabel('Temperature', fontsize=18)
    plt.ylabel('Vapor Pressure', fontsize=18)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    lgd = plt.legend(by_label.values(), by_label.keys(),fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('VLE', fontsize=20)

    # Save the plot as a PNG file
    plt.savefig('vapor_pressure.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}


process Plot_GOMC_GEMC_Production_VLE_Per_Density_VP { 
    cache 'lenient'
    fair true
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/plot_gemc/vapor_pressure", mode: 'copy', overwrite: false
    cpus 1
    debug false
    input: path(data_csv)
    output: path ("per_density.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv("$data_csv", sep='\t')
    # Create a dictionary to dynamically map each method to a specific marker, line pattern, fill pattern, and color
    method_properties = {method: (marker, linestyle, fillstyle, color) for method, marker, linestyle, fillstyle, color in zip(
        df.columns.str.split('_').str[1] + '_' + df.columns.str.split('_').str[2],
        ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h','s',],
        ['-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':','-'],  # Different line styles        
        ['full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full', 'none', 'full'],  # Alternating fill pattern
        ['red', 'red', 'blue', 'blue', 'green', 'yellow', 'orange', 'green', 'green', 'grey', 'purple'],  # Alternating colors
    )}
    method_data = {}
    # Plot for Box 0
    plt.figure(figsize=(10, 6))

    # Plot each column as a line graph with the appropriate color, marker, line pattern, and fill pattern for Box 0
    for idx, col in enumerate(df.columns):
        temperature = col.split('_')[0]
        method = col.split('_')[1] + '_' + col.split('_')[2]
        # Append temperature and density to the corresponding method in the dictionary
        if method not in method_data:
            method_data[method] = {'temperature': [], 'density': []}
        method_data[method]['temperature'].append(int(temperature))
        method_data[method]['density'].append(df[col].mean())
        if (method == "RAHBARI_DSF"):
            method = "EWALD"
            if method not in method_data:
                method_data[method] = {'temperature': [], 'density': []}
            method_data[method]['temperature'].append(int(temperature))
            method_data[method]['density'].append(df[col][0])

    # Sort temperature and density lists for each method by density
    for method, data in method_data.items():
        temp_density = list(zip(data['temperature'], data['density']))
        temp_density.sort(key=lambda x: x[1])  # Sort by density
        method_data[method]['temperature'], method_data[method]['density'] = zip(*temp_density)

    import matplotlib.pyplot as plt

    # Get unique density values
    unique_densities = sorted(set(method_data["EWALD"]["density"]))

    # Create subplots based on the number of unique density values
    num_subplots = len(unique_densities)
    fig, axs = plt.subplots(1, num_subplots, figsize=(15, 5))

    markers = ['o', 's', '^', 'D', 'v', '<', ]
    colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple'],  # Alternating colors

    # Iterate over each subplot
    # Iterate over each subplot
    for i, density in enumerate(unique_densities):
        axs[i].set_title(f'Ewald Density: {density}')
        axs[i].set_ylabel('Density')
        axs[i].tick_params(axis='x', rotation=90)  # Adjust rotation angle as needed
        # Get temperature and method values for the current density
        methods = []
        densities = []
        ew_dens = 0
        for method, data in method_data.items():
            if method != "EWALD":
                methods.append(method)
                densities.append(data['density'][i])
            else:
                ew_dens = data['density'][i]
        # Plotting bar chart for each method
        #axs[i].bar(methods, densities, color='blue')
        axs[i].scatter(methods, densities)
        axs[i].axhline(y=ew_dens, color="black", linestyle='--', linewidth=2)
    # Save the plot as a PNG file
    plt.savefig('per_density.png', bbox_inches='tight')

    # Display the plot
    plt.show()
    """
}




workflow build_NVT_system_calibrate {
    take:
    statepoint_and_solvent_xml
    jinja_channel
    main:
    build_solvent_system(statepoint_and_solvent_xml)
    write_namd_confs(build_solvent_system.out.system,jinja_channel)
    NVT_input=build_solvent_system.out.system.join(write_namd_confs.out.namd, by:0)
    NAMD_NVT_equilibration(NVT_input)
    WGC_input=build_solvent_system.out.charmm.join(NAMD_NVT_equilibration.out.restart_files, by:0)
    write_gomc_calibration_confs(WGC_input)
    GCE_input=build_solvent_system.out.system.join(NAMD_NVT_equilibration.out.restart_files,by:0).join(write_gomc_calibration_confs.out.ewald_calibration_conf,by:0)
    GOMC_Ewald_Calibration(GCE_input)
    pginput=NAMD_NVT_equilibration.out.restart_files.join(GOMC_Ewald_Calibration.out.grids, by:0)
    plot_grids(pginput)

    emit:
    //restart_files = NAMD_NVT_equilibration.out.restart_files
    //system = build_solvent_system.out.system
    //charmm = build_solvent_system.out.charmm
    convergence = plot_grids.out.convergence

    
}


workflow build_NVT_system {
    take:
    statepoint_and_solvent_xml
    jinja_channel
    main:
    build_solvent_system(statepoint_and_solvent_xml)
    write_namd_confs(build_solvent_system.out.system,jinja_channel)
    NVT_input=build_solvent_system.out.system.join(write_namd_confs.out.namd, by:0)
    NAMD_NVT_equilibration(NVT_input)
    emit:
    restart_files = NAMD_NVT_equilibration.out.restart_files
}

workflow build_GEMC_system {
    take:
    convergenceChannel
    main:
    build_two_box_system(convergenceChannel)
    write_gemc_ewald_confs(build_two_box_system.out.system)
    GOMC_GEMC_Production_Input_Channel=write_gemc_ewald_confs.out.gemc_conf
    GOMC_GEMC_Production_Replica(GOMC_GEMC_Production_Input_Channel)
    
    Extract_Density_GOMC_GEMC_Production(GOMC_GEMC_Production_Replica.out.record)
    Extract_Density_GOMC_GEMC_Production.out.analysis | collect | Collate_GOMC_GEMC_Production
    Plot_GOMC_GEMC_Production_Ewald(Collate_GOMC_GEMC_Production.out)
    //Plot_GOMC_GEMC_Production_VLE(Collate_GOMC_GEMC_Production.out)
    //Plot_GOMC_GEMC_Production_VLE_Per_Density(Collate_GOMC_GEMC_Production.out)

    //Extract_Vapor_Pressure_GOMC_GEMC_Production(GOMC_GEMC_Production_Replica.out.record)
    //Extract_Vapor_Pressure_GOMC_GEMC_Production.out.analysis | collect | Collate_GOMC_GEMC_Production_VP
    //Plot_GOMC_GEMC_Production_VLE_VP(Collate_GOMC_GEMC_Production_VP.out)
    emit:
    restart_files = GOMC_GEMC_Production_Replica.out.restart_files
}


workflow build_GEMC_system_Calibrate {
    take:
    convergenceChannel
    main:
    build_two_box_system_calibrate(convergenceChannel)
    GOMC_GEMC_Production_Input_Channel=build_two_box_system_calibrate.out.system
    GOMC_GEMC_Calibration(GOMC_GEMC_Production_Input_Channel)
    pginput = GOMC_GEMC_Calibration.out.grids
    plot_grids_two_box(pginput)
    emit:
    convergence = plot_grids_two_box.out.convergence
}


workflow build_GEMC_system_wolf {
    take:
    convergenceChannel
    main:
    methods = Channel.of( "RAHBARI_DSF","RAHBARI_DSP","WAIBEL2018_DSF","WAIBEL2018_DSP","WAIBEL2019_DSF","WAIBEL2019_DSP" )
    combinedChannel=convergenceChannel.combine(methods)
    build_two_box_system_already_calibrated(combinedChannel)
    GOMC_GEMC_Production_Input_Channel = build_two_box_system_already_calibrated.out.system
    GOMC_GEMC_Production(GOMC_GEMC_Production_Input_Channel)
    Extract_Density_GOMC_GEMC_Production(GOMC_GEMC_Production.out.record)
    Extract_Density_GOMC_GEMC_Production.out.analysis | collect | Collate_GOMC_GEMC_Production
    Plot_GOMC_GEMC_Production(Collate_GOMC_GEMC_Production.out)
    Plot_GOMC_GEMC_Production_VLE(Collate_GOMC_GEMC_Production.out)
    Plot_GOMC_GEMC_Production_VLE_Per_Density(Collate_GOMC_GEMC_Production.out)

    Extract_Vapor_Pressure_GOMC_GEMC_Production(GOMC_GEMC_Production.out.record)
    Extract_Vapor_Pressure_GOMC_GEMC_Production.out.analysis | collect | Collate_GOMC_GEMC_Production_VP
    Plot_GOMC_GEMC_Production_VLE_VP(Collate_GOMC_GEMC_Production_VP.out)
    //Plot_GOMC_GEMC_Production_VLE_VP(Collate_GOMC_GEMC_Production_VP.out)
    //Plot_GOMC_GEMC_Production_VLE_Per_Density_VP(Collate_GOMC_GEMC_Production_VP.out)

    //Plot_GEMC_Input = GOMC_GEMC_Production.out.record.collect()
    //Collate_GOMC_GEMC_Production(Plot_GEMC_Input)
}
