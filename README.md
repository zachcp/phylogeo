### [phylogeo](http://zachcp.github.com/phylogeo/)

####Geographic Mapping of Ecological Datasets

The phylogeo package is intented to allow the rapid geographic exploration of [phyloseq] (https://github.com/joey711/phyloseq) objects. The package currently only has five major: 

- `map_phyloseq()` 
- `map_network()` 
- `map_tree()`
- `map_clusters()`
- `plot_distance()`

To learn more about useage, please see the [documentation](http://zachcp.github.io/phylogeo/)

To install, phyloseq must be first be installed. phylogeo can then be installed using `devtools`:

```
install.packages("devtools")
library("devtools")
install_github("zachcp/phylogeo")
```

Maps like this one can be made with phyloseq-geo ![ExampleMap](https://dl.dropboxusercontent.com/u/1803062/examplemap.png)

