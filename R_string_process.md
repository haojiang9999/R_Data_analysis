## R string process
#### split string vector
sapply(strsplit(colnames(TM.all.matrx.filter),"_"), `[`, 3)

##### Extracting a String Between 2 Characters in R
###### By Eric Cai - The Chemical Statistician
https://chemicalstatistician.wordpress.com/2015/06/18/how-to-extract-a-string-between-2-characters-in-r-and-sas/
```
# clear all variables in workspace
rm(list=ls(all=TRUE))

# create a vector of 3 example dates
mydate = c('Jan 23/2012', 'Aug 5/2011', 'Dec 17/2011')

# getstr() is my customized function
# it extracts a string between 2 characters in a string variable
getstr = function(mystring, initial.character, final.character)
{

     # check that all 3 inputs are character variables
     if (!is.character(mystring))
     {
          stop('The parent string must be a character variable.')
     }

     if (!is.character(initial.character))
     {
          stop('The initial character must be a character variable.')
     }


     if (!is.character(final.character))
     {
          stop('The final character must be a character variable.')
     }



     # pre-allocate a vector to store the extracted strings
     snippet = rep(0, length(mystring))



     for (i in 1:length(mystring))
     {
          # extract the initial position
          initial.position = gregexpr(initial.character, mystring[i])[[1]][1] + 1

          # extract the final position
          final.position = gregexpr(final.character, mystring[i])[[1]][1] - 1

          # extract the substring between the initial and final positions, inclusively
          snippet[i] = substr(mystring[i], initial.position, final.position)
     }

     return(snippet)
}


# use the getstr() function to extract the day between the comma and the slash in "mydate"
getstr(mydate, ' ', '/')
```
