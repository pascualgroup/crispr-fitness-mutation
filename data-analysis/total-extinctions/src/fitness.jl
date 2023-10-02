function fitness(allele,centerallele,maxallele,scale,maxfitness,death,evofunction)
    if evofunction == 1
        minallele = 2*centerallele - maxallele
        scale = -1*maxfitness/((centerallele - minallele)*(centerallele - maxallele))
        return -1*(allele - minallele)*(allele - maxallele)*scale
    end
    if evofunction == 2
        return exp(-1 * scale * (allele - centerallele)^2) + death
    end
end