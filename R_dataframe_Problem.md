### How to convert entire dataframe to numeric while preserving decimals?
https://stackoverflow.com/questions/26391921/how-to-convert-entire-dataframe-to-numeric-while-preserving-decimals
```
df1[] <- lapply(df1, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(df1, class)
#         a         b 
# "numeric" "numeric" 
```
