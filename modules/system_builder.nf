#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process build_solvent_system_from_torch {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${density}", mode: 'copy', overwrite: true

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
    MC_steps=1500

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
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/input", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(temp_K), val(P_bar), val(No_mol), val(Rho_kg_per_m_cubed), val(L_m_if_cubed), path(path_to_xml)
    output:
    tuple val(Rho_kg_per_m_cubed), path("system.pdb"), path("system.psf"), path("system.inp"), emit: system
    tuple val(Rho_kg_per_m_cubed), path("statepoint.json"), emit: charmm
    path("namd_system.inp"), emit: namd_ff
    path("system_npt.conf"), emit: gomc_npt_conf
    path("ewald_calibration.conf"), emit: gomc_ewald_calibration_conf

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

    # Create a Pydantic object
    point_obj = Point(density=${Rho_kg_per_m_cubed}, temperature=${temp_K}, pressure=${P_bar},no_mol=${No_mol},box_length=${L_m_if_cubed})

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

    #MC_steps=10000
    MC_steps=1500

    gomc_control.write_gomc_control_file(charmm, conf_filename='system_nvt',  ensemble_type='NVT', RunSteps=MC_steps, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": True,
                                                            "EqSteps" : 1000,
                                                            "AdjSteps":10,
                                                            "Pressure" : float(pressure), 
                                                            "PRNG": int(0),
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "MultiParticleFreq":1.00,
                                                            "OutputName":"system_nvt"

                                                                }
                                        )
    
    restart_coor = "nvt_equil.restart.coor"
    restart_xsc = "nvt_equil.restart.xsc"
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
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.01,
                                                            "MultiParticleFreq":0.99,
                                                            "OutputName":"system_npt"

                                                                }
                                        )

    
    
    restart_coor = "system_npt_BOX_0_restart.coor"
    restart_xsc = "system_npt_BOX_0_restart.xsc"
    restart_chk = "system_npt_restart.chk"
    gomc_control.write_gomc_control_file(charmm, conf_filename='ewald_calibration',  ensemble_type='NPT', RunSteps=MC_steps, Restart=True, \
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
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.01,
                                                            "MultiParticleFreq":0.99,
                                                            "OutputName":"ewald_calibration"

                                                                }
                                        )

    file1 = open("ewald_calibration.conf", "a")
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file=restart_chk)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="1000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start="0.0",\
    end="0.5",delta="0.1")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start="10",\
    end="15",delta="0.5")
    file1.writelines(defAlphaLine)

    """
}


process write_namd_confs {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/namd_input", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path("statepoint.json")
    tuple path(path_to_minimization_template), path(path_to_nvt_template), path(path_to_npt_template)
    output:
    //tuple val(Rho_kg_per_m_cubed), path("statepoint.json"),path("system_nvt.conf"), path("system_npt.conf"), path("system.pdb"), path("system.psf"), path("system.inp"), emit: system
    tuple path("namd_min.conf"), path("namd_nvt.conf"), path("namd_npt.conf"), emit: namd

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

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.model_validate_json(json_data)

    loaded_point = load_point_from_json("statepoint.json")
    print("Loaded point")
    print(loaded_point)


    coords = {}

    minsteps=10000
    nvt_eq_steps=50000
    npt_eq_steps=100000
    #minsteps=1000
    #nvt_eq_steps=1000
    #npt_eq_steps=1000
    coords["waterModel"]="TIP3"
    coords["temp"]=loaded_point.temperature
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
    container "${params.container__namd}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/namd_npt_eq", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(statepoint),path(pdb), path(psf), path(inp), path(namd_inp)
    tuple path(namd_minimization_conf), path(namd_nvt_conf), path(namd_npt_conf)
    output:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(pdb), path(psf), path(inp), path("npt_equil.restart.xsc"), path("npt_equil.restart.coor"), emit: system
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
    container "${params.container__namd}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/namd_nvt_eq", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(pdb), path(psf), path(inp)
    path(namd_inp)
    tuple path(namd_minimization_conf), path(namd_nvt_conf), path(namd_npt_conf)
    output:
    tuple path("nvt_equil.restart.xsc"), path("nvt_equil.restart.coor"), emit: restart_files
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
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_npt_eq", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(pdb), path(psf), path(inp)
    tuple path(restart_xsc), path(restart_coor)
    path(npt_conf)
    output:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(pdb), path(psf), path(inp),\
    path("system_npt_restart.chk"),path("system_npt_BOX_0_restart.coor"), path("system_npt_BOX_0_restart.xsc"), emit: system
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
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_nvt_eq", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(statepoint),path(nvt_conf), path(npt_conf), path(pdb), path(psf), path(inp)
    output:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(npt_conf), path(pdb), path(psf), path(inp),\
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
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_npt_eq_input", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(json)
    tuple path(restart_xsc), path(restart_coor)
    output:
    tuple val(Rho_kg_per_m_cubed),path("statepoint.json"),path("system.pdb"), path("system.psf"), path("system.inp"), emit: system
    path("system_npt.conf"), emit: npt_equil_conf
    path("ewald_calibration.conf"), emit: calibration_conf
    script:
    """
    #!/usr/bin/env python
    import pickle
    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        density: float
        temperature: float
        pressure: float
        no_mol: float
        box_length: float

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


    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)

    ### Note:  The electrostatics and Ewald are turned off in the
    # GOMC control file (i.e., False) since the n-alkanes beads in the
    # trappe-ua force field have no charge (i.e., the bead charges are all zero)
    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()

    #MC_steps=10000
    MC_steps=1500

    gomc_control.write_gomc_control_file(charmm, conf_filename='system_nvt',  ensemble_type='NVT', RunSteps=MC_steps, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": True,
                                                            "EqSteps" : 1000,
                                                            "AdjSteps":10,
                                                            "Pressure" : float(pressure), 
                                                            "PRNG": int(0),
                                                            "Rcut": 12,
                                                            "RcutLow": 1,
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "MultiParticleFreq":1.00,
                                                            "OutputName":"system_nvt"

                                                                }
                                        )
    
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
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.01,
                                                            "MultiParticleFreq":0.99,
                                                            "OutputName":"system_npt"

                                                                }
                                        )

    
    
    restart_coor = "system_npt_BOX_0_restart.coor"
    restart_xsc = "system_npt_BOX_0_restart.xsc"
    restart_chk = "system_npt_restart.chk"
    gomc_control.write_gomc_control_file(charmm, conf_filename='ewald_calibration',  ensemble_type='NPT', RunSteps=MC_steps, Restart=True, \
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
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.01,
                                                            "MultiParticleFreq":0.99,
                                                            "OutputName":"ewald_calibration"

                                                                }
                                        )

    file1 = open("ewald_calibration.conf", "a")
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file=restart_chk)
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="1000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start="0.0",\
    end="0.5",delta="0.1")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start="10",\
    end="15",delta="0.5")
    file1.writelines(defAlphaLine)



    """
}

process NPT_equilibration_solvent_system {
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_npt_eq", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(npt_conf), path(pdb), path(psf), path(inp),\
    path("system_nvt_restart.chk"),path("system_nvt_BOX_0_restart.coor"), path("system_nvt_BOX_0_restart.xsc")
    output:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(npt_conf), path(pdb), path(psf), path(inp),\
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
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_ewald_calibration_input", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(json), path(charmm)
    tuple path(restart_xsc), path(restart_coor), path(restart_chk)
    output:
    tuple val(Rho_kg_per_m_cubed),path("statepoint.json"),path("system.pdb"), path("system.psf"), path("system.inp"), path("ewald_calibration.conf"), path(restart_xsc), path(restart_coor), path(restart_chk), emit: system
    tuple val(Rho_kg_per_m_cubed), path(json), path(charmm), emit: charmm
    script:
    """
    #!/usr/bin/env python
    import pickle
    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        density: float
        temperature: float
        pressure: float
        no_mol: float
        box_length: float

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


    # Build the Charmm object, which is required to write the
    # FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    # Unpickling the object
    with open("${charmm}", 'rb') as file:
        charmm = pickle.load(file)

    ### Note:  The electrostatics and Ewald are turned off in the
    # GOMC control file (i.e., False) since the n-alkanes beads in the
    # trappe-ua force field have no charge (i.e., the bead charges are all zero)
    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()

    #MC_steps=10000
    MC_steps=1500

    restart_coor = "${restart_coor}"
    restart_xsc = "${restart_xsc}"
    gomc_control.write_gomc_control_file(charmm, conf_filename='ewald_calibration',  ensemble_type='NPT', RunSteps=MC_steps, Restart=True, \
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
                                                            "RcutCoulomb_box_0": 12,
                                                            "VDWGeometricSigma" : False,
                                                            "CoordinatesFreq" : [False,1000],
                                                            "DCDFreq" : [False,100],
                                                            "RestartFreq":[True,MC_steps],
                                                            "CheckpointFreq":[True,MC_steps],
                                                            "ConsoleFreq":[True,1],
                                                            "BlockAverageFreq":[True,MC_steps],
                                                            "HistogramFreq":[False,1000],
                                                            "DisFreq":0.0,
                                                            "RotFreq":0.0,
                                                            "IntraSwapFreq":0.0,
                                                            "RegrowthFreq":0.0,
                                                            "CrankShaftFreq":0.0,
                                                            "VolFreq":0.01,
                                                            "MultiParticleFreq":0.99,
                                                            "OutputName":"ewald_calibration"

                                                                }
                                        )

    file1 = open("ewald_calibration.conf", "a")
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="Checkpoint", val="True",file="${restart_chk}")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{box}\\t{val}\\t{file}\\n".format(box="WolfCalibrationFreq", val="True",file="1000")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfAlphaRange", box="0",start="0.0",\
    end="0.5",delta="0.1")
    file1.writelines(defAlphaLine)
    defAlphaLine = "{title}\\t{box}\\t{start}\\t{end}\\t{delta}\\n".format(title="WolfCutoffCoulombRange", box="0",start="10",\
    end="15",delta="0.5")
    file1.writelines(defAlphaLine)

    """
}


process GOMC_Ewald_Calibration {
    container "${params.container__gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_ewald_calibration", mode: 'copy', overwrite: true
    cpus 8
    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(statepoint), path(pdb), path(psf), path(inp), path(npt_conf), path(restart_xsc), path(restart_coor), path(restart_chk)
    output:
    path("Wolf_Calibration_*"), emit: grids
    path("Ewald_Calibration.log"),  emit: record
    shell:
    """
    
    #!/bin/bash
    cat ${npt_conf} > local.conf
    GOMC_CPU_NPT +p${task.cpus} local.conf > Ewald_Calibration.log
    """
}

process plot_grids {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}/gomc_ewald_calibration_plots", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(Rho_kg_per_m_cubed), path(json), path(charmm)
    path(grids)
    output:
    path("grids.png"), emit: fig
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import re

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

    # Flatten the axes array to simplify indexing
    axes = axes.flatten()

    # Iterate over files and plot each DataFrame in a subplot
    for i, file_name in enumerate(file_list):
        # Extract model name
        model_name = extract_model_name(file_name)

        # Load the CSV file into a DataFrame
        df = pd.read_csv(file_name, index_col=0)

        # Plotting one line per row in the subplot
        df.plot(ax=axes[i], marker='o', linestyle='-')

        # Set plot labels and title
        axes[i].set_xlabel('Alpha')
        axes[i].set_ylabel('Relative Error (Elec Energy)')
        axes[i].set_title(f'Model: {model_name}')

    # Adjust layout
    plt.tight_layout()

    # Save the figure as 'grids.png'
    plt.savefig('grids.png')

    # Show the plot
    plt.show()
    """
}



workflow build_system {
    take:
    statepoint_and_solvent_xml
    jinja_channel
    main:
    build_solvent_system(statepoint_and_solvent_xml)
    write_namd_confs(build_solvent_system.out.charmm,jinja_channel)
    NAMD_NVT_equilibration(build_solvent_system.out.system, build_solvent_system.out.namd_ff, write_namd_confs.out.namd)
    //write_gomc_confs(build_solvent_system.out.charmm,NAMD_NVT_equilibration.out.restart_files)
    //GOMC_NPT_equilibration(write_gomc_confs.out.system, NAMD_NVT_equilibration.out.restart_files, write_gomc_confs.out.npt_equil_conf)
    //write_gomc_calibration_confs(write_gomc_confs.out.charmm,GOMC_NPT_equilibration.out.restart_files)
    //GOMC_Ewald_Calibration(write_gomc_calibration_confs.out.system)
    //GOMC_Ewald_Calibration.out.grids.view()
    //plot_grids(write_gomc_calibration_confs.out.charmm,GOMC_Ewald_Calibration.out.grids)


    emit:
    system = build_solvent_system.out.system
    
}
