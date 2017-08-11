for pop in temperate tropical
do
	gunzip ${pop}_nd/${pop}_nd.*.thetas.gz
	head -1 ${pop}_nd/${pop}_nd.aa.thetas > ${pop}_nd.thetas
	sed -i 's/#//' ${pop}_nd.thetas
	tail -qn +2 ${pop}_nd/${pop}_nd.*.thetas >> ${pop}_nd.thetas
done

mkdir window_nd