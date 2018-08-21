
new <- merge(assocs, sample_data(pruned)[,'host_genus'], by.x='sample', by.y=0, all.x=T, all.y=F)

new$rel <- new$count/new$sample_sum

pl <- ggplot(new, aes(x=sample, y=rel, fill=otu))
pl + geom_bar(color="black", stat="identity", position="stack") + facet_wrap(~ host_genus + geographic_area, scales='free_x')

