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
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), path("initial_scikit_optimize_model.pkl"), val(0), emit: scikit_optimize_model
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


process ask_points {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize/${density}/batch/${iteration}/ask", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), path(scikit_optimize_model), val(iteration)
    output:
    //tuple val(density), path(conf), path(pdb), path(psf), path(inp), path("current_scikit_optimize_model.pkl"), val(iteration), emit: scikit_optimize_model
    path("current_scikit_optimize_model.pkl"), emit: mdl
    path("*.json"), emit: json
    script:
    """
    #!/usr/bin/env python
    # for python3
    import sys
    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha: float

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
            point_obj = Point(alpha=value[0])

            # Serialize the Pydantic object to JSON
            with open(f"{count}.json", 'w') as file:
                file.write(point_obj.json())

        with open("current_scikit_optimize_model.pkl", 'wb') as f:
            pickle.dump(opt, f)
        f.close()

    """
}


process run_gomc {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/scikit_optimize/${density}/batch/${iteration}/run_gomc", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), val(iteration), path(json)
    output:
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), val(iteration), path(json)

    script:
    """
    #!/usr/bin/env python

    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        alpha: float

    # Function to load Pydantic objects from JSON file
    def load_point_from_json(file_path: str) -> Point:
        with open(file_path, 'r') as file:
            json_data = file.read()
            return Point.parse_raw(json_data)

    loaded_point = load_point_from_json("${json}")
    print("Loaded point")
    print(loaded_point)

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
    //input_csv = file(params.database_path)
    // Create a channel with the CSV file
    //csv_channel = channel.fromPath(input_csv)
    ask_points(scikit_optimize_model)
    // val(density), path(conf), path(pdb), path(psf), path(inp), path(scikit_optimize_model), val(iteration)
    density = scikit_optimize_model.map { it[0] }.view{}
    conf = scikit_optimize_model.map { it[1] }.view{}
    pdb = scikit_optimize_model.map { it[2] }.view{}
    psf = scikit_optimize_model.map { it[3] }.view{}
    inp = scikit_optimize_model.map { it[4] }.view{}
    current_iteration = scikit_optimize_model.map { it[6] }.view{}
    next_iteration = scikit_optimize_model.map { it[6]+1 }.view{}

    systems = density.merge(conf, pdb, psf, inp, current_iteration) //.view()
    systems_to_run = systems.combine(ask_points.out.json.flatten())
    run_gomc(systems_to_run)
    scikit_optimize_model = density.merge(conf, pdb, psf, inp, ask_points.out.mdl, next_iteration) //.view()

    //create_systems(ask.out.points.flatten())
    emit:
    scikit_optimize_model
    //scikit_optimize_model = ask_points.out.scikit_optimize_model
    //points = ask.out.points

}


workflow calibrate_wrapper {
    take:
    scikit_optimize_model
    main:
    calibrate(scikit_optimize_model)
    //calibrate.out.scikit_optimize_model.view()
    //calibrate.recurse(f).times(3)
    //emit:
    //points = ask.out.points

}