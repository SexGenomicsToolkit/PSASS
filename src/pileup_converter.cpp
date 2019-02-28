#include "pileup_converter.h"

PileupConverter::PileupConverter(int argc, char *argv[]) {

    if (argc != 4) {

        this->usage();
        exit(1);

    }

    if (std::string(argv[2]) == "-") {

        this->from_stdin = true;

    } else {

        this->input_file.open(std::string(argv[2]));

        if (not this->input_file.is_open()) {

            std::cerr << "Error: cannot open input file (" << std::string(argv[2]) << ")." << std::endl;
            exit(1);

        }

    }

    this->output_file.open(std::string(argv[3]));

    if (not this->output_file.is_open()) {

        std::cerr << "Error: cannot open output file (" << std::string(argv[3]) << ")." << std::endl;
        exit(1);

    }

}



void PileupConverter::usage() {

    std::cout << std::endl << "Usage: psass convert <input_file> <output_file>" << std::endl << std::endl ;
    std::cout << "Options:" << std::endl << std::endl;
    std::cout << "input-file      <string>    Either a samtools pileup output file or \"-\" for stdin" << std::endl;
    std::cout << "output-file     <string>    Output file " << std::endl;

    std::cout << std::endl;

}



void PileupConverter::process_line() {

    this->output_file << this->contig << "\t" << this->position << "\t" << this->ref_allele;
    for (uint i=0; i<6; ++i) this->output_file << "\t" << this->pool[0][i];
    for (uint i=0; i<6; ++i) this->output_file << "\t" << this->pool[1][i];
    this->output_file << "\n";

    this->field = 0;
    this->temp = "";

    for (uint i=0; i<6; ++i) {
        this->pool[0][i] = 0;
        this->pool[1][i] = 0;
    }


}



void PileupConverter::process_field() {

    switch (this->field) {

        case 0:
            this->contig = this->temp;
            break;

        case 1:
            this->position = fast_stoi(this->temp.c_str());
            break;

        case 2:
            this->ref_allele = this->temp.c_str()[0];
            break;

        case 3:
            this->depth[0] = fast_stoi(this->temp.c_str());
            break;

        case 6:
            this->depth[1] = fast_stoi(this->temp.c_str());
            break;

        default:
            break;

    }

    this->temp = "";
    ++this->field;

}



void PileupConverter::process_base(char base, bool pool) {

    if (this->read_begin) {

        this->read_begin = false;

    } else if (not this->next_indel and this->remaining_indel == 0) {

        switch (base) {

            case '.':
                ++this->pool[pool][this->bases[this->ref_allele]];
                break;

            case ',':
                ++this->pool[pool][this->bases[this->ref_allele]];
                break;

            case '*':
                ++this->pool[pool][this->bases['I']];
                break;

            case '-':
                this->next_indel = true;
                this->tmp_indel_size = "";
                ++this->pool[pool][this->bases['I']];
                break;

            case '+':
                this->next_indel = true;
                this->tmp_indel_size = "";
                ++this->pool[pool][this->bases['I']];
                break;

            case '^':
                this->read_begin = true;
                break;

            case '$':
                break;

            default:
                ++this->pool[pool][this->bases[base]];
                break;
        }

    } else if (this->next_indel) {

        if (isdigit(base)) {

            this->tmp_indel_size += base;

        } else {

            this->next_indel = false;
            this->remaining_indel = std::stoi(this->tmp_indel_size) - 1;

        }

    } else if (this->remaining_indel > 0) {

        --this->remaining_indel;

    }

}



// Read the input file and process each line
void PileupConverter::run() {

    do {

        if (not this->from_stdin) {

            this->input_file.read(this->buff, this->buff_size);
            this->k = this->input_file.gcount();

        } else {

            std::cin.read(this->buff, this->buff_size);
            this->k = std::cin.gcount();

        }

        for (uint i=0; i<this->k; ++i) {

            switch (this->buff[i]) {

                case '\r':
                    break;

                case '\n':
                    this->process_line();
                    break;

                case '\t':
                    this->process_field();
                    break;

                default:

                    switch (field) {

                        case 4:
                            this->process_base(this->buff[i], 0);
                            break;

                        case 7:
                            this->process_base(this->buff[i], 1);
                            break;

                        case 5:
                            break;

                        case 8:
                            break;

                        default:
                            this->temp += this->buff[i];
                            break;
                    }

                    break;

            }
        }

    } while ((not this->from_stdin and this->input_file) or (this->from_stdin and std::cin));

}
