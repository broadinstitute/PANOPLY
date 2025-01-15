library(KEGGgraph)
library(igraph)



p = "hsa00620" # pyruvate, arbitrarily chosen pathway
# g = pathways_qs$graph.list[[p]] # pull pathway as is from pathway_qs, as a framework for what I want to make
# pdf('/opt/input/tmp.pdf'); with(environment(MetaboAnalystR::RerenderMetPAGraph), plotGraph(g)); dev.off()

# retrieve KGML file from KEGG database
fn = paste0(p,'.kgml')
KEGGgraph::retrieveKGML(p, 'hsa', destfile = fn)

# # parse file directly to graph
# kg = parseKGML2Graph(fn,expandGenes=TRUE)

# parse to pathway and THEN to graph
mapkpathway = parseKGML(fn)
kg = KEGGpathway2Graph(mapkpathway,
                       genesOnly = F,
                       expandGenes=F)
kgp = KEGGpathway2reactionGraph(mapkpathway)
pdf('/opt/input/plotKEGGreaction.pdf'); plot(kgp); dev.off()

# # parse node metadata from pathway
# map_list = sapply(mapkpathway@nodes, function(n) {
#   list(entryID = n@entryID,
#        id = n@name,
#        type = n@type,
#        link = n@link,
#        reaction = n@reaction,
#        #map = n@map,
#        name = n@graphics@name,
#        x = n@graphics@x,
#        y = n@graphics@y,
#        shape = n@graphics@type,
#        width = n@graphics@width,
#        height = n@graphics@height,
#        fgcolor = n@graphics@fgcolor,
#        bgcolor = n@graphics@bgcolor)
# }) # parse info out of S4 KEGG Node objects


#### Convert from Pathway to Graph ####
nchar_limit_map = 25
nchar_limit = 6
# pull node data out of mapkpathway@nodes
vertex_df = data.frame(name = sapply(mapkpathway@nodes, function(n) {n@entryID}),
                       id = sapply(mapkpathway@nodes, function(n) {n@name %>% paste(collapse="")}),
                       type = sapply(mapkpathway@nodes, function(n) {n@type}),
                       link = sapply(mapkpathway@nodes, function(n) {n@link}),
                       reaction = sapply(mapkpathway@nodes, function(n) {n@reaction}),
                       name_hr = sapply(mapkpathway@nodes, function(n) {n@graphics@name}),
                       x = sapply(mapkpathway@nodes, function(n) {n@graphics@x}),
                       y = -sapply(mapkpathway@nodes, function(n) {n@graphics@y}),
                       shape = sapply(mapkpathway@nodes, function(n) {n@graphics@type}),
                       # size = sapply(mapkpathway@nodes, function(n) {n@graphics@width}),
                       # size2 = sapply(mapkpathway@nodes, function(n) {n@graphics@height}),
                       frame.color = sapply(mapkpathway@nodes, function(n) {n@graphics@fgcolor}),
                       color = sapply(mapkpathway@nodes, function(n) {n@graphics@bgcolor})) %>%
  # add vertex formatting (https://igraph.org/r/doc/plot.common.html#:~:text=The%20list%20of%20parameters)
  dplyr::mutate(shape = ifelse(shape=="roundrectangle", "rectangle", shape),
                color = ifelse(type=="map", "grey90", color),
                frame.color = ifelse(type=="map", "grey90", frame.color),
                label = ifelse(type=="map", stringr::str_trunc(name_hr, nchar_limit_map),
                               stringr::str_trunc(name_hr, nchar_limit)), # 
                # label = stringr::str_trunc(name_hr, nchar_limit),
                # label.distance = ifelse(),
                size = ifelse(type=="compound", 3, 15*(nchar(label)/8)),
                label.dist = ifelse(type=="compound", .75, 0),
                size2 = 8)
# vertex_df[c('x','y')] = norm_coords(as.matrix(vertex_df[c('x','y')]))

# pull edge data out of mapkpathway@edges
edge_df = data.frame(e1 = sapply(mapkpathway@edges, function(e) {e@entry1ID}),
                     e2 = sapply(mapkpathway@edges, function(e) {e@entry2ID}),
                     type = sapply(mapkpathway@edges, function(e) {e@type}),
                     st.n = sapply(mapkpathway@edges, function(e) {e@subtype$subtype@name}),
                     st.v = sapply(mapkpathway@edges, function(e) {e@subtype$subtype@value})) %>%
  # add edge formatting
  dplyr::mutate(label = st.v,
                lty = 1, # by default, solid line
                lty = ifelse(e1 %in% dplyr::filter(vertex_df, type=="map")$name, 2, lty), # make lines from maps dashed
                lty = ifelse(e2 %in% dplyr::filter(vertex_df, type=="map")$name, 2, lty), # make lines to maps dashed
                curved = FALSE)
# add edge to intermediary compound
edge_full_df = data.frame('e1' = c(edge_df$e1, edge_df$st.v), # add intermediate edge to compound
                          'e2' = c(edge_df$st.v, edge_df$e2)) %>% # add edge from compound to gene
  dplyr::mutate(lty = 1, # by default, solid line
                lty = ifelse(e1 %in% dplyr::filter(vertex_df, type=="map")$name, 2, lty), # make lines from maps dashed
                lty = ifelse(e2 %in% dplyr::filter(vertex_df, type=="map")$name, 2, lty)) # make lines to maps dashed
# v_missing = setdiff(rownames(vertex_df), unique(c(edge_df$e1, edge_df$e2)))
# mg = graph_from_data_frame(edge_df, vertices = vertex_df)
mg = graph_from_data_frame(edge_full_df, vertices = vertex_df) %>%
  simplify(edge.attr.comb = "first")
# plot graph
pdf(glue('/opt/input/{p}.pdf'), width = 20, height=20)
plot.igraph(mg)
dev.off()
# pdf(glue('/opt/input/{p}_nicely.pdf'), width = 16, height=16)
# plot.igraph(mg, layout=layout_nicely(mg))
# dev.off()
pdf(glue('/opt/input/{p}_noMaps.pdf'), width = 16, height=16)
mg_noMap = delete_vertices(mg, vertex_df[vertex_df$type=="map", 'name'])
plot.igraph(mg_noMap, edge.label = E(mg_noMap)$st.v)
dev.off()
pdf(glue('/opt/input/{p}_simplify.pdf'), width = 16, height=16)
mg_simple = delete.vertices(mg, degree(mg)==0) # delete vertices without edges
plot.igraph(mg_simple, edge.label = E(mg_simple)$st.v)
dev.off()

# mg = make_empty_graph(n = dim(map_df)[1], directed = TRUE) %>%
#   set_vertex_attr("name", value = rownames(map_df)) %>%
#   set_vertex_attr("name", value = kg@nodes) %>%
  
  
#### Exploring Features missing from MetaboAnalyst Pathways ####
# determine what nodes are excluded from "final" graph
included_nodes = map_list[( names(map_list) %in% names(V(g)) )]
excluded_nodes = map_list[!( names(map_list) %in% names(V(g)) )]

# what type are the in/excluded nodes
lapply(included_nodes, function(x) {x$type}) %>% unique() # only genes/compounds are included
lapply(excluded_nodes, function(x) {x$type}) %>% unique() # some genes and compounds ARE excluded though...?

excluded_entries = excluded_nodes[unlist(lapply(excluded_nodes, function(x) {x$type %in% c('gene', 'compound')}))]
sapply(included_nodes, function(x) {x$name}) 
sapply(excluded_entries, function(x) {x$name}) 


#### Direct igraph Conversion ####
ig = igraph.from.graphNEL(kg) # missing a ton of info that would be needed for nice plotting
pdf('/opt/input/tmp4.pdf'); plot.igraph(ig); dev.off()

ig.s = delete.vertices(simplify(ig), degree(ig)==0) # delete vertices without edges


# pdf('/opt/input/tmp.pdf'); with(environment(MetaboAnalystR::RerenderMetPAGraph), plotGraph(g)); dev.off()