#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Using recursion
nextflow.preview.recursion=true

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { build_system } from './modules/system_builder'
include { train_model } from './modules/model_builder'
include { predict_model } from './modules/model_builder'
include { initialize_scikit_optimize_model } from './modules/scikit_optimize'
include { calibrate_wrapper } from './modules/scikit_optimize'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run [-profile local/docker/singularity/slurm] . [--database_path /Path/To/CSV --num_points num_points_coex_curve]

# FreeSolv
nextflow run . -profile docker

# User defined database
nextflow run . -profile docker --database_path databases/spce.csv

Required Arguments:

  Input Data:
    --database_path       /Path/To/CSV.
Optional Arguments:
    --output_folder       Folder for output files (default $projectDir/results).

    """.stripIndent()
}


// Main workflow
workflow {


log.info """\
=============================================================================
         output_folder          : ${params.output_folder}
         database_path          : ${params.database_path}
         path_to_xml            : ${params.path_to_xml}
         batch_size             : ${params.batch_size}
         density_lb             : ${params.density_lb}
         density_ub             : ${params.density_ub}
         alpha_lb               : ${params.alpha_lb}
         alpha_ub               : ${params.alpha_ub}
         torch_model            : ${params.torch_model}
         torch_scalers          : ${params.torch_scalers}


         """
         .stripIndent()


    // Supported charge types : 'am1bcc', 'am1-mulliken', 'gasteiger', 'resp'
    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.database_path == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }
    if ( params.database_path ){
        // Define the input CSV file
        input_csv = file(params.database_path)
        //densities = Channel.fromList( ['1', '10', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'] )
        densities = Channel.fromList( ['1', '10'] )

        // Create a channel with the CSV file
        csv_channel = channel.fromPath(input_csv)
        solventData = Channel.fromPath( params.database_path ).splitCsv(header: true,limit: 2,quote:'"').map { 
            row -> [row.temp_K, row.P_bar, row.No_mol, row.Rho_kg_per_m_cubed, row.L_m_if_cubed]
        }
        //vapor_systems = build_solvents(vapor_points.combine(path_to_xml))
        //path_to_database = Channel.fromPath( params.database_path )
        //train_model(csv_channel)
        //model_density_tuple = train_model.out.model_scalers_tuple.combine(densities)
        //predicted_points = predict_model(model_density_tuple)
        solvent_xml = file(params.path_to_xml)
        solvent_xml_channel = channel.fromPath(solvent_xml)
        system_input = solventData.combine(solvent_xml_channel)
        system_input.view()
        build_system(system_input)
        skopt_model = initialize_scikit_optimize_model(build_system.out.system)
        skopt_model.view()
        calibrate_wrapper(skopt_model)
    } else {
        helpMessage()
    }
}