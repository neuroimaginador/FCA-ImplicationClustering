# planets dataset - Clustering of Implications
## Setup

First, we load the needed libraries


```r
# Suppress excess info
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)
library("tidyverse")
library("cluster")
library("fcaR")
```

Load other functions needed to compute implication
distances and perform clustering.


```r
source("./distances_registry.R")
source("./implication_distance.R")
source("./cluster.R")
source("./one_hot.R")
source("./count_attributes.R")
```

Auxiliary functions


```r
imps2str <- function(imp) {
  sapply(seq(imp$cardinality()), 
         function(j) 
           fcaR:::.implication_to_string(
             imp[j]$get_LHS_matrix(), 
             imp[j]$get_RHS_matrix(), 
             imp$get_attributes()))
}
```


 
## The planets dataset

We'll use the planets dataset, which is included
in the fcaR package, to show how to cluster 
implications.

First, load and print the dataset.



```r
library(fcaR)
fc_planets <- FormalContext$new(planets)

fcaR:::.print_binary(planets) %>% 
  kableExtra::kable(align = "c")
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> small </th>
   <th style="text-align:center;"> medium </th>
   <th style="text-align:center;"> large </th>
   <th style="text-align:center;"> near </th>
   <th style="text-align:center;"> far </th>
   <th style="text-align:center;"> moon </th>
   <th style="text-align:center;"> no_moon </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Mercury </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Venus </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Earth </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mars </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Jupiter </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Saturn </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Uranus </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Neptune </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pluto </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;"> X </td>
   <td style="text-align:center;">  </td>
  </tr>
</tbody>
</table>

Let us find the Duquenne-Guigues basis of
implications and make a copy in variable "imps"



```r
fc_planets$find_implications()
imps <- fc_planets$implications$clone()
imps
```

```
## Implication set with 10 implications.
## Rule 1: {no_moon} -> {small, near}
## Rule 2: {far} -> {moon}
## Rule 3: {near} -> {small}
## Rule 4: {large} -> {far, moon}
## Rule 5: {medium} -> {far, moon}
## Rule 6: {medium, large, far, moon} -> {small,
##   near, no_moon}
## Rule 7: {small, near, moon, no_moon} ->
##   {medium, large, far}
## Rule 8: {small, near, far, moon} -> {medium,
##   large, no_moon}
## Rule 9: {small, large, far, moon} -> {medium,
##   near, no_moon}
## Rule 10: {small, medium, far, moon} -> {large,
##   near, no_moon}
```


We can compute the dissimilarity between each pair
of implications, as
$\text{diss}(P\to Q, R\to T) := d(P, R) + d(P^+, R^+)$
where $d$ is the Manhattan distance.



```r
Distances <- impl_dist(imps, method = "Manhattan")
Dimps <- Distances$left + Distances$closure
# DD is the dissiilarity matrix
DD <- as.matrix(Dimps)
kableExtra::kable(DD, align = "c", row.names = TRUE)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> 1 </th>
   <th style="text-align:center;"> 2 </th>
   <th style="text-align:center;"> 3 </th>
   <th style="text-align:center;"> 4 </th>
   <th style="text-align:center;"> 5 </th>
   <th style="text-align:center;"> 6 </th>
   <th style="text-align:center;"> 7 </th>
   <th style="text-align:center;"> 8 </th>
   <th style="text-align:center;"> 9 </th>
   <th style="text-align:center;"> 10 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
</tbody>
</table>


We can use a clustering algorithm to find $k = 2$
clusters


```r
cluster <- cluster_implications(Dimps, k = 2, "pam")

# The central implications
imps[cluster$representatives]
```

```
## Implication set with 2 implications.
## Rule 1: {far} -> {moon}
## Rule 2: {small, medium, far, moon} -> {large,
##   near, no_moon}
```


The implications in the first cluster are



```r
imps[cluster$clustering == 1]
```

```
## Implication set with 5 implications.
## Rule 1: {no_moon} -> {small, near}
## Rule 2: {far} -> {moon}
## Rule 3: {near} -> {small}
## Rule 4: {large} -> {far, moon}
## Rule 5: {medium} -> {far, moon}
```


and in the second one:



```r
imps[cluster$clustering == 2]
```

```
## Implication set with 5 implications.
## Rule 1: {medium, large, far, moon} -> {small,
##   near, no_moon}
## Rule 2: {small, near, moon, no_moon} ->
##   {medium, large, far}
## Rule 3: {small, near, far, moon} -> {medium,
##   large, no_moon}
## Rule 4: {small, large, far, moon} -> {medium,
##   near, no_moon}
## Rule 5: {small, medium, far, moon} -> {large,
##   near, no_moon}
```

Note that this cluster is composed by all implications
in which all attributes appear. Let us remove them
to cluster the other implications properly:


```r
imps <- imps[imps$support() > 0]
print(imps)
```

```
## Implication set with 5 implications.
## Rule 1: {no_moon} -> {small, near}
## Rule 2: {far} -> {moon}
## Rule 3: {near} -> {small}
## Rule 4: {large} -> {far, moon}
## Rule 5: {medium} -> {far, moon}
```


Let us compute the dissimilarity matrix



```r
imps
```

```
## Implication set with 5 implications.
## Rule 1: {no_moon} -> {small, near}
## Rule 2: {far} -> {moon}
## Rule 3: {near} -> {small}
## Rule 4: {large} -> {far, moon}
## Rule 5: {medium} -> {far, moon}
```

```r
Distances <- impl_dist(imps, method = "Manhattan")
Dimps <- Distances$left + Distances$closure
DD <- as.matrix(Dimps)
kableExtra::kable(DD, align = "c", row.names = TRUE)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> 1 </th>
   <th style="text-align:center;"> 2 </th>
   <th style="text-align:center;"> 3 </th>
   <th style="text-align:center;"> 4 </th>
   <th style="text-align:center;"> 5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
</tbody>
</table>


And cluster into two groups of implications



```r
cluster <- cluster_implications(Dimps, k = 2, "pam")

# The central implications are
imps[cluster$representatives]
```

```
## Implication set with 2 implications.
## Rule 1: {near} -> {small}
## Rule 2: {far} -> {moon}
```


Now, the implications in the first cluster are:



```r
imps[cluster$clustering == 1]
```

```
## Implication set with 2 implications.
## Rule 1: {no_moon} -> {small, near}
## Rule 2: {near} -> {small}
```


and the ones in the second one are:



```r
imps[cluster$clustering == 2]
```

```
## Implication set with 3 implications.
## Rule 1: {far} -> {moon}
## Rule 2: {large} -> {far, moon}
## Rule 3: {medium} -> {far, moon}
```


This way, we can tag the clusters according to the
closed sets that arise from the central implications



```r
for (i in seq(k)) {
  
  cat(" - Cluster", i, ": ")
  imp_med <- imps[cluster$representatives[i]]$clone()
  M <- fcaR:::.union(imp_med$get_LHS_matrix(),
                     imp_med$get_RHS_matrix())
  S <- SparseSet$new(attributes = fc_planets$attributes,
                     M = M)
  print(S)
  cat("\n")
  
}
```

 - Cluster 1 : {small, near}
 - Cluster 2 : {far, moon}


This shows that the dataset has two main clusters:

1) Planets near the Sun and therefore small,

2) Planets far from the Sun and therefore with satellites. 

__Graphical representation of the implication space__

We can use multidimensional scaling to plot the 
implications and their clusters:


```r
imp_mds <- MASS::isoMDS(d = Dimps, k = 2, trace = FALSE) 
df <- imp_mds$point %>% as.data.frame() %>% 
  as_tibble() %>%  
  set_names(c("x", "y")) %>% 
  mutate(cluster = cluster$clustering,
         imp = imps2str(imps))

library(ggrepel)
ggplot(df, aes(x = x, y = y, color = cluster, label = imp)) + 
  geom_label_repel() +
  geom_point() +
  theme_void() + 
  theme(legend.position = "none",
        axis.line = element_line(size = 0))
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)


__What happens if we consider only the left-hand side
of the implications when computing the dissimilarity
matrix?__



```r
Distances <- impl_dist(imps, method = "Manhattan")
Dimps <- Distances$left
DD <- as.matrix(Dimps)
rownames(DD) <- seq(imps$cardinality())
kableExtra::kable(DD, align = "c", row.names = TRUE)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> 1 </th>
   <th style="text-align:center;"> 2 </th>
   <th style="text-align:center;"> 3 </th>
   <th style="text-align:center;"> 4 </th>
   <th style="text-align:center;"> 5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
</tbody>
</table>


We can use the PAM clustering algorithm to find 
2 clusters:



```r
cluster <- cluster_implications(Dimps, k = 2, "pam")

# The central implications in this case are:
imps[cluster$representatives]
```

```
## Implication set with 2 implications.
## Rule 1: {large} -> {far, moon}
## Rule 2: {medium} -> {far, moon}
```


The implications in the first cluster:



```r
imps[cluster$clustering == 1]
```

```
## Implication set with 4 implications.
## Rule 1: {no_moon} -> {small, near}
## Rule 2: {far} -> {moon}
## Rule 3: {near} -> {small}
## Rule 4: {large} -> {far, moon}
```


and in the second one:



```r
imps[cluster$clustering == 2]
```

```
## Implication set with 1 implications.
## Rule 1: {medium} -> {far, moon}
```

This clustering is almost random, due to the uninformative 
dissimilarity matrix used.
