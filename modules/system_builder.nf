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

    MC_steps=10000
    gomc_control.write_gomc_control_file(charmm, conf_filename='system',  ensemble_type='NVT', RunSteps=MC_steps, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": False,
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
    publishDir "${params.output_folder}/systems/density_${Rho_kg_per_m_cubed}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(temp_K), val(P_bar), val(No_mol), val(Rho_kg_per_m_cubed), val(L_m_if_cubed), path(path_to_xml)
    output:
    tuple val(Rho_kg_per_m_cubed), path("statepoint.json"),path("system_nvt.conf"), path("system_npt.conf"), path("system.pdb"), path("system.psf"), path("system.inp"), emit: system

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
        file.write(point_obj.json())

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
                                        n_compounds=int(${No_mol}+int(float(${No_mol})*0.1)),
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

    MC_steps=10000

    gomc_control.write_gomc_control_file(charmm, conf_filename='system_nvt',  ensemble_type='NVT', RunSteps=MC_steps, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": False,
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

    gomc_control.write_gomc_control_file(charmm, conf_filename='system_npt',  ensemble_type='NPT', RunSteps=MC_steps, Restart=True, \
                                        check_input_files_exist=False, Temperature=float(temperature) * u.Kelvin, ExpertMode=True,\
                                        Coordinates_box_0="system.pdb",Structure_box_0="system.psf",binCoordinates_box_0="system_nvt_BOX_0_restart.coor",
                                        extendedSystem_box_0="system_nvt_BOX_0_restart.xsc",
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": False,
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




    """
}


workflow build_system {
    take:
    statepoint_and_solvent_xml
    main:
    build_solvent_system(statepoint_and_solvent_xml)
    
    emit:
    system = build_solvent_system.out.system
    
}
