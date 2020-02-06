#include "pileup_converter.h"

PileupConverter::PileupConverter(Parameters& parameters) {

    /* Constructor for PileupConverter.
     * Handles arguments parsing for the 'convert' subcommand.
     */

    if (parameters.input_file_path == "-") {
        this->from_stdin = true;
    } else {
        this->input_file.open(parameters.input_file_path);
        this->input_file_path = parameters.input_file_path;
    }

    if (parameters.output_file_path != "") {
        this->to_stdout = false;
        this->ofile.open(parameters.output_file_path);
        this->output_file_path = parameters.output_file_path;
    }

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
            i == 5 ? this->obuff[this->j] = '\t' : this->obuff[this->j] = ',';
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
            i == 5 ? this->obuff[this->j] = '\t' : this->obuff[this->j] = ',';
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
            i == 5 ? this->obuff[this->j] = '\t' : this->obuff[this->j] = ',';
            ++this->j;

        } else if (this->pool[pool_number][i] > 9) {

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i] / 10);
            this->pool[pool_number][i] = this->pool[pool_number][i] % 10;
            ++this->j;
            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            i == 5 ? this->obuff[this->j] = '\t' : this->obuff[this->j] = ',';
            ++this->j;

        } else {

            this->obuff[this->j] = '0' + static_cast<char>(this->pool[pool_number][i]);
            ++this->j;
            i == 5 ? this->obuff[this->j] = '\t' : this->obuff[this->j] = ',';
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

    --this->j;
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
    this->contig = "";

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
            if (this->contig != this->current_contig or this->current_contig == "") {
                if (this->current_contig == "") this->current_contig = this->contig;
                std::string header = "region=" + this->current_contig + "\t len=NA";
                if (not this->to_stdout) {
                    this->ofile.write(this->current_contig.c_str(), static_cast<uint>(this->current_contig.size())) ;  // Write contig to output file
                } else {
                    fwrite(this->current_contig.c_str(), static_cast<uint>(this->current_contig.size()), 1, stdout);   // Write contig to sdout
                }
                this->current_contig = this->contig;
            }
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

    log("Pileup converter started");

    this->from_stdin ? log("Reading input from stdin") : log("Reading input from <" + this->input_file_path + ">");
    this->to_stdout ? log("Outputting to stdout") : log("Outputting to <" + this->output_file_path + ">");

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
                            this->contig += this->buff[i];  // Field 0 --> scaffold name
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

    log("Pileup converter ended successfully");
}
