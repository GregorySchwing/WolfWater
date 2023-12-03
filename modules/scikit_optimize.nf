#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process initialize_model {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize/${density}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(density), path(conf), path(pdb), path(psf), path(inp)
    output:
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), path("initial_scikit_optimize_model.pkl"), emit: scikit_optimize_model
    path("log.txt")
    script:
    """
    #!/usr/bin/env python
    # for python3
    import sys
    with open("log.txt", 'w') as sys.stdout:
        from skopt import Optimizer
        import pickle
        import numpy as np
        np.int = np.int64
        opt = Optimizer([(${params.alpha_lb},${params.alpha_ub})],
                    "GP", acq_func="EI",
                    acq_optimizer="sampling",
                    initial_point_generator="lhs",
                    n_initial_points=${params.batch_size},
                    random_state=123)

        with open('initial_scikit_optimize_model.pkl', 'wb') as f:
            pickle.dump(opt, f)
        f.close()
        print("Model initialized with domain=",[(${params.alpha_lb},${params.alpha_ub})])
        print("GP")
        print("acq_func=EI")
        print("acq_optimizer=sampling")
        print("initial_point_generator=lhs"),
        print("n_initial_points/batch size=",${params.batch_size})


    """
}


process ask {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize/batch/${task.index}/ask", mode: 'copy', overwrite: true

    debug true
    input:
    path(scikit_optimize_model)
    output:
    path("current_scikit_optimize_model.pkl"), emit: scikit_optimize_model
    path("*.json"), emit: points
    script:
    """
    #!/usr/bin/env python
    # for python3
    import sys
    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha: float
        density: float

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
            point_obj = Point(alpha=value[0], density=value[1])

            # Serialize the Pydantic object to JSON
            with open(f"{count}.json", 'w') as file:
                file.write(point_obj.json())

        with open("current_scikit_optimize_model.pkl", 'wb') as f:
            pickle.dump(opt, f)
        f.close()

    """
}


process create_systems {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/scikit_optimize/batch/${task.index}/systems", mode: 'copy', overwrite: true

    debug true
    input:
    path(json)
    output:

    script:
    """
    #!/usr/bin/env python

    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha: float
        density: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.parse_raw(json_data)

    loaded_point = load_point_from_json("${json}")
    print("Loaded point")
    print(loaded_point)
    quit()
    print("#**********************")
    print("Started: GOMC Charmm Object")
    print("#**********************")

    total_molecules_liquid = 500
    box_0_box_size_ang = 25

    forcefield_files = '../../SPCE_GMSO.xml'
    molecule_A = mb.load('../../SPCE.mol2')
    molecule_A.name = 'WAT'

    molecule_type_list = [molecule_A]
    mol_fraction_molecule_A = 1.0
    molecule_mol_fraction_list = [mol_fraction_molecule_A]
    fixed_bonds_angles_list = [molecule_A.name]
    residues_list = [molecule_A.name]

    #bead_to_atom_name_dict = {'_CH3': 'C', '_CH2': 'C', '_CH': 'C', '_HC': 'C'}

    print('total_molecules_liquid = ' + str(total_molecules_liquid))

    print('Running: liquid phase box packing')
    box_liq = mb.fill_box(compound=molecule_type_list,
                          n_compounds=total_molecules_liquid,
                          box=[box_0_box_size_ang/10, box_0_box_size_ang/10, box_0_box_size_ang/10]
                          )
    # this uses the UFF force field, not the user force field.  Causes trouble in simulations of rigid molecules
    #box_liq.energy_minimize(forcefield=forcefield_files,
    #                        steps=10 ** 5
    #                        )
    print('Completed: liquid phase box packing')

    print('molecule_mol_fraction_list =  ' + str(molecule_mol_fraction_list))

    print('Running: GOMC FF file, and the psf and pdb files')
    gomc_charmm = mf_charmm.Charmm(
        box_liq,
        mosdef_structure_box_0_name_str,
        ff_filename=gomc_ff_filename_str,
        forcefield_selection=forcefield_files,
        residues=residues_list,
        gomc_fix_bonds_angles=fixed_bonds_angles_list,
    )

    if write_files == True:
        gomc_charmm.write_inp()

        gomc_charmm.write_psf()

        gomc_charmm.write_pdb()

    print('Completed: Writing  GOMC FF file, and the psf and pdb files')
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
    input_csv = file(params.database_path)
    // Create a channel with the CSV file
    csv_channel = channel.fromPath(input_csv)
    ask(scikit_optimize_model)
    create_systems(ask.out.points.flatten())
    emit:
    scikit_optimize_model = ask.out.scikit_optimize_model
    //points = ask.out.points

}
