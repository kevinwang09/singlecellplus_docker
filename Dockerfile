# DO NOT EDIT FILES CALLED 'Dockerfile'; they are automatically
# generated. Edit 'Dockerfile.in' and generate the 'Dockerfile'
# with the 'rake' command.

# The suggested name for this image is: bioconductor/release_base.

FROM bioconductor/release_core2

MAINTAINER kevin.wang@sydney.edu.au

ADD install.R /home/


RUN wget https://storage.googleapis.com/scp_data/SingleCellPlus_data.zip -P /home/rstudio/
RUN cd /home/rstudio/ && unzip SingleCellPlus_data.zip
RUN mkdir /home/rstudio/data
RUN cp -r /home/rstudio/SingleCellPlus_data/* /home/rstudio/data
RUN rm -rf /home/rstudio/SingleCellPlus_data.zip
RUN rm -rf /home/rstudio/SingleCellPlus_data
RUN git clone https://github.com/SydneyBioX/SingleCellPlus
RUN cp ./SingleCellPlus/qc.Rmd /home/rstudio/
RUN cp ./SingleCellPlus/scMerge.Rmd /home/rstudio/
RUN cp ./SingleCellPlus/downstream.Rmd /home/rstudio/
RUN rm -rf /home/rstudio/SingleCellPlus
RUN ls /home/rstudio/


# RUN R -f /home/install.R
# ADD test.R /home/rstudio
# RUN ls /home/rstudio/