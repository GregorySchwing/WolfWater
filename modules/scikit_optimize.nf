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
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), val(iteration), path("*.json"), emit: systems
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
    publishDir "${params.output_folder}/scikit_optimize/${density}/batch/${iteration}/run_gomc/${json.baseName}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(density), path(conf), path(pdb), path(psf), path(inp), val(iteration), path(json)
    output:
    tuple val(density), path("system.conf"), path(pdb), path(psf), path(inp), val(iteration), path(json)

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
    # Specify the file paths
    input_file_path = "${conf}"  # Replace with your actual input file path
    output_file_path = "system.conf"  # Replace with your desired output file path
    print(input_file_path)
    print(output_file_path)
    # Open input file for reading and output file for writing
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        # Read the input file line by line
        for line in input_file:
            # Write each line to the output file
            output_file.write(line)
        output_file.write("FUCK!")

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
    ask_points.out.systems.transpose().view()
    run_gomc(ask_points.out.systems.transpose())

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