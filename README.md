# MCB432_MA-2
#Db4 script for R studio- 3D Scatter plot graph
# Read aligned sequences (MSA) into R as a DNAbin object
db4<- read.FASTA("SeqDatabase_db4_unique.edit.filtered.aln.fasta", type = "DNA")
#Identify the structure of the file
class(db4)
typeof(db4)
length(db4)
#Generate DNA distance with dist.dna() JC69 is a default evolution/substitution model
dist_db4<- dist.dna(db4,model="JC69")
#check the structure of the distance file
str(dist_db4)
#Use phangorn package to generate a neighbor-joining tree
db4.nj<- nj(dist_db4)
plot(db4.nj,cex=0.5); add.scale.bar(0, -2)
# Use nodelabels() to label nodes
dev.new(height=12,width=16);plot(db4.nj,cex=0.5: "Neighbor-JoiningTree"); add.scale.bar(0, -2); nodelabels(db4.nj$node)
# Identify the outgroup and use grep to match any tip labels containing "Methanobacterium"(outgroup)
outgroup_tips <- grep("Methanobacterium", db4.nj$tip.label)
# Print the matching tip labels to confirm
print(db4.nj$tip.label[outgroup_tips])
# Root the tree using the identified outgroup tips
db4.rooted <- root(db4.nj, outgroup=outgroup_tips, resolve.root=TRUE)
# Finding subtree, this command list all subtrees.
 subtrees(db4.rooted, wait=FALSE)
# Use distance matrix for ordination and scatter plots
db4.pcoa<- pcoa(dist_db4)
str(db4.pcoa)
#check the vectors (1:2- 2D)
(db4.PCoA<- db4.pcoa$vectors[,1:2])
# Create a matrix for Accession, Genus and Species and evaluate the strucutre(syntax) of the matrix
mat = matrix(nrow = 144, ncol = 3)
colnames(mat)<- c("Acc", "Gen", "sp")
for (k in 1:144){mat[k,]<- strsplit(rownames(db4.PCoA)," ")[[k]]}
head(mat)
#Plot the PCoA
plot(db4.PCoA, main="PCoA of db4"); text(db4.PCoA,label=mat[,1]
# Create a matrix for Accession, Genus and Species (cluster)
km<- kmeans(db4.PCoA,5);
# Use kcol for clusters, plot with dots and text (2d scatter plot)
plot(db4.PCoA, xlab="PCo1", ylab="PCo2", main="PCoA of db4", col=km$cluster)
layout(matrix(1:2,1,2)); km<- kmeans(db4.PCoA,5,nstart=25)
plot(db4.PCoA, main=paste("PCoA of db4, Kmeans k=5"), col=km$cluster,xlab="PCo1",
ylab="PCo2", pch=19);
plot(db4.PCoA, main=paste("PCoA of db4, Kmeans k=5"), xlab="PCo1",ylab="PCo2",type="n");
text(db4.PCoA,label=paste(mat[,2],mat[,3], sep=" "),cex=0.5, col=km$cluster)
#Establishing three principalcoordinates for the 3d-plot
# Check the structure of db4.PCoA
str(db4.PCoA)
# Check the dimensions (rows = samples, columns = PCoA axes)
dim(db4.PCoA)
# Perform PCoA with 3 dimensions (for 3D plotting)
db4.PCoA_3D <- cmdscale(dist_db4, k=3)
# Check the first few rows of the PCoA result
head(db4.PCoA_3D)
#Install plotly (3D scatter plot)
Create a 3D scatter plot with plotly (without text labels)
fig <- plot_ly(
x = db4.PCoA_3D[, 1],   # First principal coordinate (PCo1)
y = db4.PCoA_3D[, 2],   # Second principal coordinate (PCo2)
z = db4.PCoA_3D[, 3],   # Third principal coordinate (PCo3)
type = 'scatter3d', 
mode = 'markers',  # Just markers, no text
marker = list(size = 5, 
color = km$cluster,  # Color by K-means cluster colorscale = 'Viridis', opacity = 0.8))
# Add axis labels and a plot title
fig <- fig %>% layout(title = "3D PCoA of db4 with K-means Clustering (k=5)",scene = list(xaxis = list(title = 'PCo1'),yaxis = list(title = 'PCo2'),zaxis = list(title = 'PCo3')))
# Show the plot
fig
# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(fig, "PCoA_Kmeans_3D.html")

