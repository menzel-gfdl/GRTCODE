#ifndef PARSE_CSV_H_
#define PARSE_CSV_H_


/*Parse a csv file.*/
int parse_csv(char const * const filepath,
              int * const num_lines,
              int * const num_cols,
              char *** out);


#endif
