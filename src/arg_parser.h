#pragma once
#include <iostream>
#include <stdio.h>
#include "CLI11/CLI11.hpp"
#include "parameters.h"

// Failure message function for CLI parser
inline std::string failure_message(const CLI::App* parser, const CLI::Error& error) {

    std::string message = "";

    if (error.what() == std::string("A subcommand is required")) {
        message = "\nSubcommand error: missing or invalid subcommand\n\n" + parser->help();
    } else if (error.get_exit_code() == 106) {  // 106 corresponds to wrong argument type
        message = "\nArgument error: " + std::string(error.what()) + "\n\n" + parser->help();
    } else {
        message = "\nError: " + std::string(error.what()) + "\n\n" + parser->help();
    }

    return message;
}


// Formatter for CLI
class CustomFormatter : public CLI::Formatter {

    public:

        uint column_widths[3] {0, 0, 0};  // Will be used to store the maximum width of each column : flags, type, description
        uint border_width = 4;  // Space between two columns

        // Formatter for an Option line, overrides the same function from CLI::Formatter
        virtual std::string make_option(const CLI::Option* opt, bool is_positional) const {

            std::string option = "", name = "", type = "", description = "", default_value = "", required = "REQUIRED";

            // Generate option name, if positional -> just the name, if not positional -> <short_flag, long_flag>
            name = opt->get_name();
            type = opt->get_type_name();
            description = opt->get_description();
            default_value = opt->get_default_str();

            // Generate the help string for this option, adding the right number of spaces after each column based on column_widths
            option = name + std::string(border_width + column_widths[0] - name.size(), ' ');
            option += type + std::string(border_width + column_widths[1] - type.size(), ' ');
            option += description + std::string(border_width + column_widths[2] - description.size(), ' ');
            if (opt->get_required()) default_value = required;
            if (default_value != "") option += "[" + default_value + "]";
            option += "\n";

            return option;
        }

        virtual std::string make_description(const CLI::App *app) const {

            return "";
        }

        void set_column_widths(CLI::App& parser) {
            std::string tmp = "";
            for (auto opt: parser.get_options()) {
                opt->get_positional() ? tmp = opt->get_name() : tmp = "--" + opt->get_lnames()[0];
                if (tmp.size() > this->column_widths[0]) this->column_widths[0] = static_cast<uint>(tmp.size());
                tmp = opt->get_type_name();
                if (tmp.size() > this->column_widths[1]) this->column_widths[1] = static_cast<uint>(tmp.size());
                tmp = opt->get_description();
                if (tmp.size() > this->column_widths[2]) this->column_widths[2] = static_cast<uint>(tmp.size());
            }
        }

};


// Argument parsing main function
inline Parameters parse_args(int& argc, char** argv) {

    CLI::App parser {""};  // Parser instance from CLI App parser
    Parameters parameters;

    std::shared_ptr<CustomFormatter> formatter(new CustomFormatter);

    // Main parser options
    parser.formatter(formatter);  // Set custom help format defined above
    parser.require_subcommand();  // Check that there is a subcommand
    parser.failure_message(failure_message);  // Formatting for error message

    CLI::App* analyze = parser.add_subcommand("analyze", "Compute metrics from a sync file from psass convert or from popoolation2.");
    CLI::App* convert = parser.add_subcommand("convert", "Convert a pileup file from samtools to a synchronized pool file.");

    // Options for 'analyze'
    analyze->add_option("INPUT_FILE", parameters.input_file_path, "Path to a sync file generated by psass convert or popoolation2")->required()->check(CLI::ExistingFile);
    analyze->add_option("OUTPUT_FILE", parameters.output_file_path, "Path to an output file for sliding window metrics")->required();

    analyze->add_option("--pool1", parameters.pool1_id, "Name of the first pool (order in the pileup file)", true)->group("Input/Output");
    analyze->add_option("--pool2", parameters.pool2_id, "Name of the second pool (order in the pileup file)", true)->group("Input/Output");
    analyze->add_option("--gff-file", parameters.gff_file_path, "Path to a GFF file for gene-specific output", true)->group("Input/Output");
    analyze->add_flag("--popoolation", parameters.popoolation_format, "If set, assumes the input file was generated with popoolation2")->group("Input/Output");
    analyze->add_option("--snps-file", parameters.snp_pos_file_path, "Output sex-biased SNPs to this file", true)->group("Input/Output");
    analyze->add_option("--fst-file", parameters.fst_pos_file_path, "Output high FST positions to this file", true)->group("Input/Output");

    analyze->add_option("--min-depth", parameters.min_depth, "Minimum depth to include a site in the analyses", true)->group("Analysis");
    analyze->add_option("--window-size", parameters.window_size, "Size of the sliding window (in bp)", true)->group("Analysis");
    analyze->add_option("--output-resolution", parameters.output_resolution, "Output resolution for sliding window metrics (in bp)", true)->group("Analysis");
    analyze->add_option("--freq-het", parameters.freq_het, "Frequency of a sex-linked SNP in the heterogametic sex", true)->group("Analysis");
    analyze->add_option("--range-het", parameters.range_het, "Range of frequency for a sex-linked SNP in the heterogametic sex", true)->group("Analysis");
    analyze->add_option("--freq-hom", parameters.freq_hom, "Frequency of a sex-linked SNP in the homogametic sex", true)->group("Analysis");
    analyze->add_option("--range-hom", parameters.range_hom, "Range of frequency for a sex-linked SNP in the homogametic sex", true)->group("Analysis");
    analyze->add_option("--min-fst", parameters.min_fst, "Minimum FST to output a site in the FST positions file", true)->group("Analysis");
    analyze->add_flag("--group-snps", parameters.group_snps, "If set, group consecutive snps to count them as a single polymorphism")->group("Analysis");

    // Options for 'convert'
    convert->add_option("INPUT", parameters.input_file_path, "Either a path to a samtools pileup output file or \"-\" for stdin")->required();
    //    option->check(CLI::ExistingFile);  // NEED TO ADD A CUSTOM CHECK IF VALUE IS NOT -

    convert->add_option("--output-file", parameters.output_file_path, "Write to an output file instead of stdout");

    // The parser throws an exception upon failure and implements an exit() method which output an error message and returns the exit code.
    try {

        parser.parse(argc, argv);

    } catch (const CLI::ParseError &e) {

        if (parser.get_subcommands().size() > 0) {

            std::string tmp = "";
            for (auto opt: parser.get_subcommands()[0]->get_options()) {
                tmp = opt->get_name();
                if (tmp.size() > formatter->column_widths[0]) formatter->column_widths[0] = static_cast<uint>(tmp.size());
                tmp = opt->get_type_name();
                if (tmp.size() > formatter->column_widths[1]) formatter->column_widths[1] = static_cast<uint>(tmp.size());
                tmp = opt->get_description();
                if (tmp.size() > formatter->column_widths[2]) formatter->column_widths[2] = static_cast<uint>(tmp.size());
            }

        } else {
            formatter->column_widths[0] = 10;
            formatter->column_widths[1] = 0;
            formatter->column_widths[2] = 50;
        }

        exit(parser.exit(e));

    }

    // Set some parameter values after parsing
    parameters.min_het = parameters.freq_het - parameters.range_het;
    parameters.max_het = parameters.freq_het + parameters.range_het;
    parameters.min_hom = parameters.freq_hom- parameters.range_hom;

    // Get subcommand name
    CLI::App* subcommand = parser.get_subcommands()[0];
    parameters.command = subcommand->get_name();

    return parameters;
}
