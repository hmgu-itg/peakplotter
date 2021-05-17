Bootstrap: docker
From: ubuntu:18.04

%environment
	TZ=Europe/Berlin
	PATH=$PATH:/usr/local/bin:/opt/locuszoom/bin:/opt/plink
    
    LC_ALL=C.UTF-8
    LANG=C.UTF-8

%post
    apt update
	ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
    
    apt install -y tabix moreutils git


# locuszoom    
    cd /opt
	wget https://statgen.sph.umich.edu/locuszoom/download/locuszoom_1.4_srconly.tgz
	tar -zxf locuszoom_1.4_srconly.tgz
	rm locuszoom_1.4_srconly.tgz
	
	mkdir locuszoom/conf
	touch locuszoom/conf/m2zfast.conf
	cat <<-EOF >locuszoom/conf/m2zfast.conf
	METAL2ZOOM_PATH="bin/locuszoom.R"
	NEWFUGUE_PATH=""
	PLINK_PATH="plink"
	RSCRIPT_PATH="Rscript"
	TABIX_PATH="tabix"
	SQLITE_DB={"b38":""}
	DEFAULT_BUILD="b38"
	DEFAULT_POP="EUR"
	DEFAULT_SOURCE="b38"
	GWAS_CATS={"b38":{}}
	LD_DB={"b38":{}}
	EOF


# Bedtools
    apt-get install -y zlib1g zlib1g-dev firefox python-dev emacs

    cd /opt
    git clone https://github.com/arq5x/bedtools2.git
	cd bedtools2
	# gwava needs a version without "sam header" error messages
	git checkout tags/v2.27.1
	make && make install

	pip2 install pybedtools scipy pandas numpy scikit-learn==0.14.1


# PLINK
	cd /opt
	mkdir plink
	cd plink
	wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
	unzip plink_linux_x86_64_20201019.zip
	rm plink_linux_x86_64_20201019.zip