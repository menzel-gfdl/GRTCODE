#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "utils.h"

#define MAXCHARSPERLINE 1024
#define MAXCHARSPERTOKEN 32


#define offset(col,line,nlines) (col*nlines + line)


/*Copy an input token into an input array.*/
static int copy_token(char * const token,
                      int const col_index,
                      int const line_index,
                      int const num_lines,
                      char **vals)
{
    int s = strlen(token);
    if (s <= 0 || s >= MAXCHARSPERTOKEN)
    {
        /*The input token has an invalid size.*/
        fatal(VALUE_ERR,
              "the size (%d) of the token %s on line %d must be in the"
                  " range [1,%d].",
              s,
              token,
              line_index+1,
              MAXCHARSPERTOKEN);
    }
    if (token[s-1] == '\n')
    {
        /*Strip off the trailing new line character.*/
        token[s-1] = '\0';
    }
    strncpy(vals[offset(col_index,line_index,(num_lines))],
            token,
            (size_t)MAXCHARSPERTOKEN);
    return SUCCESS;
}


/*Parse a csv file, assuming that the first line in the file contains
  headers for each of the columns.*/
int parse_csv(char const * const filepath,
              int * const num_lines,
              int * const num_cols,
              int const ignore_headers,
              char *** out)
{
    /*Check inputs.*/
    not_null(filepath);
    not_null(num_lines);
    not_null(num_cols);
    not_null(out);

    /*Open the file.*/
    FILE *f = fopen(filepath,"r");
    if (f == NULL)
    {
        fatal(IO_ERR,
              "failed to open csv file %s.",
              filepath);
    }

    /*Count the number of lines/columns on a line.*/
    char line[MAXCHARSPERLINE];
    *num_lines = 0;
    *num_cols = -1;
    while (fgets(line,MAXCHARSPERLINE,f) != NULL)
    {
        (*num_lines)++;
        char *c = line;
        if (*c == '\n')
        {
            fatal(VALUE_ERR,
                  "line %d in file %s is blank.",
                  *num_lines,
                  filepath);
        }
        int n = 1;
        int num_chars = 0;
        while (*c != '\n')
        {
            if (*c == ',')
            {
                n++;
            }
            c++;
            num_chars++;
            if (num_chars > MAXCHARSPERLINE)
            {
                fatal(VALUE_ERR,
                      "the number of characters (>=%d) on line %d of file %s"
                          " exceeds the maximum allowed (%d).",
                      num_chars,
                      *num_lines,
                      filepath,
                      MAXCHARSPERLINE);
            }
        }
        if (*num_cols < 0)
        {
            *num_cols = n;
        }
        else if (n != *num_cols)
        {
            fatal(VALUE_ERR,
                  "the number of columns (%d) on line %d of file"
                      " %s differs from the number of columns (%d) on"
                      " the other lines.",
                  n,
                  *num_lines,
                  filepath,
                  *num_cols);
        }
    }
    if (*num_lines == 0)
    {
        fatal(VALUE_ERR,
              "the file %s is empty.",
              filepath);
    }

    /*Allocate arrays to hold the values that will be read in.*/
    rewind(f);
    char **vals;
    if (ignore_headers)
    {
        (*num_lines)--;
        if (*num_lines == 0)
        {
            fatal(VALUE_ERR,
                  "the file %s only contains headers, no data.",
                  filepath);
        }
        fgets(line,MAXCHARSPERLINE,f);
    }
    int num_vals = (*num_cols)*(*num_lines);
    check(malloc_ptr((void **)(&vals),
                     sizeof(*vals)*num_vals));
    int i;
    for (i=0;i<num_vals;++i)
    {
        check(malloc_ptr((void **)(&(vals[i])),
                         sizeof(*(vals[i]))*MAXCHARSPERTOKEN));
        snprintf(vals[i],
                 MAXCHARSPERTOKEN,
                 "%c",
                 '\0');
    }

    /*Read in the data.*/
    int line_index = 0;
    while (fgets(line,MAXCHARSPERLINE,f) != NULL)
    {
        int col_index = 0;
        char *token = strtok(line,",");
        check(copy_token(token,
                         col_index,
                         line_index,
                         *num_lines,
                         vals));
        while (1)
        {
            token = strtok(NULL,",");
            if (token == NULL)
            {
                break;
            }
            col_index++;
            check(copy_token(token,
                             col_index,
                             line_index,
                             *num_lines,
                             vals));
        }
        line_index++;
    }
    *out = vals;

    /*Close the file.*/
    if (0 != fclose(f))
    {
        fatal(IO_ERR,
              "failed to close csv file %s.",
              filepath);
    }
    return SUCCESS;
}
