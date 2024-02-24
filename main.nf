#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Using recursion
// nextflow.preview.recursion=true

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { build_NVT_system } from './modules/system_builder'
include { build_GEMC_system } from './modules/system_builder'
include { build_GEMC_system_Calibrate } from './modules/system_builder'
include { build_GEMC_system_wolf } from './modules/system_builder'
include { build_GEMC_system_wolf_inside_vle } from './modules/system_builder'

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
         debugging              : ${params.debugging}


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
        solventData = Channel.fromPath( params.database_path ).splitCsv(header: true,limit: 14,quote:'"').map { 
            row -> [row.temp_K, row.P_bar, row.No_mol, row.Rho_kg_per_m_cubed, row.L_m_if_cubed, row.RcutCoulomb]
        }
        solvent_xml = file(params.path_to_xml)
        min_jinja_template = file(params.path_to_minimization_template)
        nvt_jinja_template = file(params.path_to_nvt_template)
        npt_jinja_template = file(params.path_to_npt_template)

        solvent_xml_channel = Channel.fromPath( solvent_xml )
        jinja_channel = Channel.fromPath( [file(params.path_to_minimization_template),\
        file(params.path_to_nvt_template), file(params.path_to_npt_template)] ).collect()
        system_input = solventData.combine(solvent_xml_channel)     
        build_NVT_system(system_input,jinja_channel)        
        convergenceChannel = build_NVT_system.out.restart_files.groupTuple(by:0,size:2,remainder:false)
        convergenceChannelFlattened = convergenceChannel.map { tuple ->
            def temperature = tuple[0]
            def densities = tuple[1]
            def statepointPaths = tuple[2]
            def xscPaths = tuple[3]
            def coorPaths = tuple[4]
            // Customize this part based on your specific requirements
            return [temperature, densities[0], densities[1], statepointPaths[0],statepointPaths[1], \
            xscPaths[0], xscPaths[1], coorPaths[0], coorPaths[1]]
        }
        gemc_system_input = convergenceChannelFlattened.combine(solvent_xml_channel)   
        build_GEMC_system(gemc_system_input)       
        tempAndDensity = convergenceChannel.map { tuple ->
            def temperature = tuple[0]
            def densities = tuple[1]
            def statepointPaths = tuple[2]
            // Customize this part based on your specific requirements
            return [temperature, densities[0], densities[1], statepointPaths[0],statepointPaths[1]]
        }
        gemc_calibration_input = tempAndDensity.join(build_GEMC_system.out.restart_files).combine(solvent_xml_channel)
        build_GEMC_system_Calibrate(gemc_calibration_input)
        gemc_wolf_production_input = gemc_system_input.join(build_GEMC_system_Calibrate.out.convergence)
        build_GEMC_system_wolf_inside_vle(gemc_wolf_production_input,build_GEMC_system.out.ewald_density_data,build_GEMC_system.out.ewald_vapor_pressure_data,build_GEMC_system.out.ewald_vol_data)
        return
        build_GEMC_system_wolf(gemc_wolf_production_input,build_GEMC_system.out.ewald_density_data,build_GEMC_system.out.ewald_vapor_pressure_data,build_GEMC_system.out.ewald_vol_data)
        return

        skopt_model = initialize_scikit_optimize_model(build_NVT_system.out.system)
        skopt_model.view()
        calibrate_wrapper(skopt_model)
    } else {
        helpMessage()
    }
}