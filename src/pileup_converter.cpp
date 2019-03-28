#include "pileup_converter.h"
#include <string.h>
#include <stdio.h>

PileupConverter::PileupConverter(int argc, char *argv[]) {

    /* Constructor for PileupConverter.
     * Handles arguments parsing for the 'convert' subcommand.
     */

    if (argc != 3 and argc != 5) {  // 3 arguments if no output file, 5 arguments if output file

        this->usage();
        exit(1);

    }

    if (std::string(argv[2]) == "-") {  // Detect if input should be read from stdin

        this->from_stdin = true;  // Boolean indicating reading from stdin

    } else {

        this->input_file.open(std::string(argv[2]));  // Try to open input file if specified

        if (not this->input_file.is_open()) {

            std::cerr << "Error: cannot open input file (" << std::string(argv[2]) << ")." << std::endl;
            exit(1);

        }
    }

    if (argc == 5) {

        this->ofile.open(argv[4]);

        if (not this->ofile.is_open()) {  // Try to open output file if specified

            std::cerr << "Error: cannot open output file (" << std::string(argv[4]) << ")." << std::endl;
            exit(1);

        }

    } else {

        this->to_stdout = true;  // Boolean indicating writing to sdout

    }

}



void PileupConverter::usage() {

    /* Simple usage output function for PileupConverter.
     */

    std::cout << std::endl << "Usage: psass convert <input_file> [ -o <output_file> ]" << std::endl << std::endl ;
    std::cout << "Options:" << std::endl << std::endl;
    std::cout << "input-file            <string>    Either a samtools pileup output file or \"-\" for stdin" << std::endl;
    std::cout << "-o <output_file>      <string>    Write output to <output_file> instead of stdout" << std::endl;
    std::cout << std::endl;

}



void PileupConverter::add_counts_to_buffer(uint8_t pool_number) {

    /* Convert integers for nucleotide counts in a pool to char and add them to the output buffer.
     * Because the values will be output in plain text, each char can only store a single digit.
     */

    for (int i=0; i<6; ++i) {

        // Code is a bit repetitive but this way only necessary calculations are performed.
        if (this->pool[pool_number][i] > 9999) {  // Max int value is 65535

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 10000);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 10000;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 1000);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 1000;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 100);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 100;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 10);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 10;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            this->obuff[this->j] = '\t';
            ++this->j;

        } else if (this->pool[pool_number][i] > 999) {

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 1000);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 1000;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 100);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 100;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 10);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 10;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            this->obuff[this->j] = '\t';
            ++this->j;

        } else if (this->pool[pool_number][i] > 99) {

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 100);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 100;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 10);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 10;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            this->obuff[this->j] = '\t';
            ++this->j;

        } else if (this->pool[pool_number][i] > 9) {

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 10);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 10;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            this->obuff[this->j] = '\t';
            ++this->j;

        } else {

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            this->obuff[this->j] = '\t';
            ++this->j;

        }
    }
}



void PileupConverter::process_line() {

    /* Process a single line from the pileup file and outputs in converted format.
     * Convert each number from the counter to a series of digits and store them in output buffer.
     */

    this->add_counts_to_buffer(0);  // Add nucleotide counts to output buffer for first pool
    this->add_counts_to_buffer(1);  // Add nucleotide counts to output buffer for second pool

    --j;
    this->obuff[this->j] = '\n';  // Add return character at the end of the line
    ++this->j;

    if (not this->to_stdout) {
        this->ofile.write(this->obuff, this->j);  // Write buffer to output file
    } else {
        fwrite(this->obuff, j, 1, stdout);  // Write buffer to sdout
    }

    // Reset variables
    this->field = 0;
    this->temp = "";

    for (uint i=0; i<6; ++i) {
        this->pool[0][i] = 0;
        this->pool[1][i] = 0;
    }

}



void PileupConverter::process_field() {

    /* Process a single field from the pileup output.
     * Add '\t' to the output buffer for relevant fields.
     */

    switch (this->field) {

        case 0:
            this->obuff[this->j] = '\t';  // Add tab after contig field
            ++this->j;
            break;

        case 1:
            this->obuff[this->j] = '\t';  // Add tab after position field
            ++this->j;
            break;

        case 2:
            this->obuff[this->j] = '\t';  // Add tab after reference base field
            ++this->j;
            break;

        default:
            break;

    }

    ++this->field;
}



void PileupConverter::process_base(char base, bool pool) {

    /* Process a single base from pileup field 4 or field 6 (nucleotides).
     * - ATGC/atgc increment respective counters
     * - ^ indicates beginning of a read and next base is mapping quality (PHRED), they are both skipped
     * - * indicates an indel that was described previously, increment indel counter
     * - + and - indicate the start of an indel with format '[+/-][length][sequence]'. Indels of any length are counted as 1.
     * Indel length L is recovered from sequence and the following L bases (corresponding to indel sequence) are skipped.
     * - $ indicates end of a read and is skipped.
     * - Any other character will trigger a warning.
     */

    if (base == '.' or base == ',') base = this->ref_allele;  // Change base to reference allele if necessary (easier processing later)

    if (this->read_begin) {  // Skip character if it's a mapping quality score.

        this->read_begin = false;

    } else if (not this->next_indel and this->remaining_indel == 0) {

        // Using a full switch here is 1x faster than using an unordered_map
        switch (base) {

            case '^':  // Character '^' marks the beginning of a read and is followed by a character indicating mapping quality. This character is skipped.
                this->read_begin = true;
                break;

            case 'A':
                ++this->pool[pool][0];
                break;

            case 'a':
                ++this->pool[pool][0];
                break;

            case 'T':
                ++this->pool[pool][1];
                break;

            case 't':
                ++this->pool[pool][1];
                break;

            case 'C':
                ++this->pool[pool][2];
                break;

            case 'c':
                ++this->pool[pool][2];
                break;

            case 'G':
                ++this->pool[pool][3];
                break;

            case 'g':
                ++this->pool[pool][3];
                break;

            case 'N':
                ++this->pool[pool][4];
                break;

            case 'n':
                ++this->pool[pool][4];
                break;

            case '*':  // * represent indels already specified before
                ++this->pool[pool][5];
                break;

            case '-':
                this->next_indel = true;  // Indels are handled separately
                this->tmp_indel_size = "";
                ++this->pool[pool][5];
                break;

            case '+':
                this->next_indel = true;  // Indels are handled separately
                this->tmp_indel_size = "";
                ++this->pool[pool][5];
                break;

            case '$':
                break;

            default:
                std::cerr << "Warning: unknown character <" << base << "> in pileup file" << std::endl;
                break;
        }

    } else if (this->next_indel) {  // Start of indel string

        if (isdigit(base)) {  // Check if character is a digit and add it to indel size string buffer if it is

            this->tmp_indel_size += base;

        } else {

            this->next_indel = false;  // If not a digit, character is the start of indel sequence string.
            this->remaining_indel = static_cast<uint>(std::stoi(this->tmp_indel_size) - 1);  // Initialize indel remaining bases variable based on size string buffer

        }

    } else if (this->remaining_indel > 0) {

        --this->remaining_indel;  // Decrement indel remaining bases

    }

}




void PileupConverter::run() {

    /* Main function from PileupConverter.
     * Read the input stream into a buffer and iterate over the buffer.
     * Reading input streams this way is much faster than using std::getline()
     */

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
                    this->j=0;
                    break;

                case '\t':
                    this->process_field();
                    break;

                default:  // non-special characters

                    switch (field) {

                        case 0:
                            this->obuff[j] = this->buff[i];  // Field 0 --> scaffold name
                            ++this->j;
                            break;

                        case 1:
                            this->obuff[j] = this->buff[i];  // Field 1 --> position on scaffold
                            ++this->j;
                            break;

                        case 2:
                            this->obuff[j] = this->buff[i];  // Field 2 --> reference allele
                            ++this->j;
                            this->ref_allele = this->buff[i];
                            break;

                        case 4:
                            this->process_base(this->buff[i], 0);  // Field 4 --> nucleotides in first pool
                            break;

                        case 7:
                            this->process_base(this->buff[i], 1);  // Field 7 --> nucleotides in second pool
                            break;

                        default:
                            break;
                    }

                    break;

            }
        }

    } while ((not this->from_stdin and this->input_file) or (this->from_stdin and std::cin));

}
