#ifndef PARSE_CSV_H_
#define PARSE_CSV_H_


/** @brief Parse a csv file, assuming that the first line in the file contains
           headers for each of the columns.
    @return RS_SUCCESS or an error code.*/
int parse_csv(char const * const filepath, /**< csv file.*/
              int * const num_lines, /**< Number of lines in the file.*/
              int * const num_cols, /**< Number of columns in the file.*/
              int const ignore_headers, /**< Flag to ignore the first line in the file.*/
              char *** out /**< Array of values.*/
             );


#endif
