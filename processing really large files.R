getwd()
setwd("~/Google Drive/post trx reg data/datafiles_screen1_miseq/from_server/")

one = read.csv("pre_one.csv")
two = read.csv("pre_two.csv")
three = read.csv("pre_three.csv")
four = read.csv("pre_four.csv")
five = read.csv("pre_five.csv")
six = read.csv("pre_six.csv")
seven = read.csv("pre_seven.csv")
eight = read.csv("pre_eight.csv")
nine = read.csv("pre_nine.csv")
ten = read.csv("pre_ten.csv")
eleven = read.csv("pre_eleven.csv")
twelve = read.csv("pre_twelve.csv")
thirteen = read.csv("pre_thirteen.csv")
fourteen = read.csv("pre_fourteen.csv")
fifteen = read.csv("pre_fifteen.csv")
sixteen = read.csv("pre_sixteen.csv")
seventeen = read.csv("pre_seventeen.csv")
eighteen = read.csv("pre_eighteen.csv")

head(eighteen)

one_pre = rbind(one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen)

table(one_pre$joint_frame)

barplot(table(one_pre$joint_frame))

write.csv(one_pre, "screen1 no recombination.csv")

pdf(sprintf("aa_start_and_end fract of cds %s.pdf", "pre_recomb_1"), useDingbats = FALSE)
hist(one_pre$dist_aa_start, main=sprintf("aa_start fract of orf %s", "pre_recomb_1"))
hist(one_pre$dist_aa_end, main=sprintf("aa_end fract of orf %s", "pre_recomb_1"))
dev.off()

#out of curiousity, see length of read vs count
pdf(sprintf("frag count v readlength %s.pdf", "pre_recomb_1"), useDingbats = FALSE)
plot(log(one_pre$frag_count, base = 10), one_pre$read_length,
     main=sprintf("frag count v length %s", "pre_recomb_1"))
dev.off()


head(one_pre)

getwd()
xref = read.delim("../SGD_features.tab", header=FALSE, quote="")
head(xref)
head(match(one_pre$gene_name, xref$V4))
head(one_pre)
median(matrix(table(one_pre$gene_name)))
table(matrix(table(one_pre$gene_name)))

write.csv()

post = read.csv("../post screen1.csv")
head(post)
table(post$gene_name)
matrix(table(post$gene_name))
table(matrix(table(post$gene_name)))
median(matrix(table(post$gene_name)))

up = read.csv("../up screen1.csv")
head(up)
table(up$gene_name)
matrix(table(up$gene_name))
median(matrix(table(up$gene_name)))
table(matrix(table(up$gene_name)))

down = read.csv("../down screen1.csv")
head(down)
table(down$gene_name)
matrix(table(down$gene_name))
median(matrix(table(down$gene_name)))

gene_count_down = matrix(table(down$gene_name))
gene_count_down
hist(gene_count_down)
table(gene_count_down)
