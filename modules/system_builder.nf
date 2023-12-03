#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process build_solvent_system {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/systems/temperature_${temperature}_density_${density}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(temperature), val(density), val(dens_uncertainty), path(path_to_xml)
    output:
    tuple val(temperature), val(density), path("system.*"), emit: system

    script:
    """
    #!/usr/bin/env python
    # GOMC Example for the NVT Ensemble using MoSDeF [1, 2, 5-10, 13-17]

    # Import the required packages and specify the force field (FF),
    # box dimensions, density, and mol ratios [1, 2, 5-10, 13-17].

    import mbuild as mb
    import unyt as u
    import mosdef_gomc.formats.gmso_charmm_writer as mf_charmm
    import mosdef_gomc.formats.gmso_gomc_conf_writer as gomc_control
    from   mbuild.lib.molecules.water import WaterSPC
    import foyer

    forcefield_file_water = "${path_to_xml}"

    liquid_box_length_Ang = 33

    liquid_box_density_kg_per_m_cubed = ${density}

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
                            ff_filename="NVT_water",
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
    gomc_control.write_gomc_control_file(charmm, conf_filename='in_NVT.conf',  ensemble_type='NVT', RunSteps=MC_steps, Temperature=float("${temperature}") * u.Kelvin, ExpertMode=True,\
                                        input_variables_dict={"ElectroStatic": True,
                                                            "Ewald": False,
                                                            "PRNG": int(0),
                                                            "Rcut": 12 * u.angstrom,
                                                            "RcutLow": 1 * u.angstrom,
                                                            "RcutCoulomb_box_0": 12 * u.angstrom,
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

workflow build_solvents {
    take:
    solv_temp_xml
    main:
    build_solvent_system(solv_temp_xml)
    
    emit:
    system = build_solvent_system.out.system
    
}
