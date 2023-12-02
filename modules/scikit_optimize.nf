#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process initialize_model {
    container "${params.container__scikit_optimize}"
    publishDir "${params.output_folder}/scikit_optimize", mode: 'copy', overwrite: true

    debug true
    input:
    output:
    path("initial_scikit_optimize_model.pkl"), emit: scikit_optimize_model
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
        opt = Optimizer([(${params.alpha_lb},${params.alpha_ub}),(${params.density_lb},${params.density_ub})],
                    "GP", acq_func="EI",
                    acq_optimizer="sampling",
                    initial_point_generator="lhs",
                    n_initial_points=${params.batch_size},
                    random_state=123)

        with open('initial_scikit_optimize_model.pkl', 'wb') as f:
            pickle.dump(opt, f)
        f.close()
        print("Model initialized with domain=",[(${params.density_lb},${params.density_ub}),(${params.alpha_lb},${params.alpha_ub})])
        print("GP")
        print("acq_func=EI")
        print("acq_optimizer=sampling")
        print("initial_point_generator=lhs"),
        print("n_initial_points/batch size=",${params.batch_size})


    """
}

workflow initialize_scikit_optimize_model {
    take:
    main:
    initialize_model()
    emit:
    scikit_optimize_model = initialize_model.out.scikit_optimize_model
    
}
