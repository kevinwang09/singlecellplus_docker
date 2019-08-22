# DO NOT EDIT FILES CALLED 'Dockerfile'; they are automatically
# generated. Edit 'Dockerfile.in' and generate the 'Dockerfile'
# with the 'rake' command.

# The suggested name for this image is: bioconductor/release_base.

FROM bioconductor/release_core2

MAINTAINER kevin.wang@sydney.edu.au

ADD install.R /home/


RUN wget http://shiny.maths.usyd.edu.au/singlecellplus/SingleCellPlus_zip.zip -P /home/rstudio/
RUN cd /home/rstudio/ && unzip SingleCellPlus_zip.zip
RUN cp -r /home/rstudio/SingleCellPlus_zip/* /home/rstudio/
RUN rm -rf /home/rstudio/SingleCellPlus_zip.zip
RUN rm -rf /home/rstudio/SingleCellPlus_zip
RUN git clone https://github.com/SydneyBioX/SingleCellPlus
RUN cp ./SingleCellPlus/qc.Rmd /home/rstudio/
RUN cp ./SingleCellPlus/scMerge.Rmd /home/rstudio/
RUN cp ./SingleCellPlus/downstream.Rmd /home/rstudio/
RUN rm -rf /home/rstudio/SingleCellPlus
RUN ls /home/rstudio/

RUN R -f /home/install.R
ADD test.R /home/rstudio
RUN cd /home/rstudio/
RUN R -f /home/rstudio/test.R