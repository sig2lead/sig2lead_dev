FROM rocker/shiny:3.5.1 

RUN apt-get update && apt-get install -y sudo libcurl4-gnutls-dev libxml2-dev libssl-dev git; 

RUN rm -r /srv/shiny-server; 

RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny; 

RUN sed -i 's/3838/3838 0.0.0.0/' /etc/shiny-server/shiny-server.conf 

RUN R -e "install.packages(pkgs=c('tidyverse','shinydashboard'), repos='https://cran.rstudio.com/')" 
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
#RUN R -e "devtools::install_github('andrewsali/shinycssloaders')"
RUN R -e "install.packages('lubridate', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('glue', repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install()"
RUN R -e "BiocManager::install(c('ChemmineOB', 'ChemmineR'))"
RUN R -e "install.packages(c('shiny', 'shinyjs','ggplot2','plotly','httr','jsonlite','DT','visNetwork','plyr','dendextend','colorspace','ggforce','rlist','scatterpie','ggrepel','bazar','XML','RCurl','bitops','scrapeR','igraph','circlize','enrichR','readr','dplyr','gplots'))"

# RUN git clone https://github.com/sig2lead/sig2lead_beta.git

# RUN cd /sig2lead_beta/Sig2Lead_beta

#copy the app to the image
# COPY ./Sig2Lead_Backend.Rproj /srv/shiny-server/
# COPY ./ui.R /srv/shiny-server/
# COPY ./sample2.txt /srv/shiny-server/
# COPY ./LINCSCompounds.csv /srv/shiny-server/
# COPY ./lincs_fps.RData /srv/shiny-server/
# COPY ./minSim_apfp_RObjects.Rdata /srv/shiny-server/
# COPY ./ui.R /srv/shiny-server/
# COPY ./server.R /srv/shiny-server/
# COPY ./lib /srv/shiny-server/
# COPY ./similarity_search /srv/shiny-server/

#COPY /sig2lead_beta/Sig2Lead_beta/ /srv/shiny-server/ 

COPY /app /srv/shiny-server/ 

RUN chmod -R +rx /srv/shiny-server/ 

EXPOSE 3838 

CMD ["/usr/bin/shiny-server.sh"]


  

# install R packages required 
# (change it dependeing on the packages you need)
# RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
# RUN R -e "devtools::install_github('andrewsali/shinycssloaders')"
# RUN R -e "install.packages('lubridate', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('glue', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('DT', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"

# RUN R -e "install.packages('BiocManager')"
# RUN R -e "BiocManager::install()"
# RUN R -e "BiocManager::install(c('ChemmineOB', 'ChemmineR'))"
# RUN R -e "install.packages(c('shiny', 'shinyjs','ggplot2','plotly','httr','jsonlite','DT','visNetwork','plyr','dendextend','colorspace','ggforce','rlist','scatterpie','ggrepel','bazar','XML','RCurl','bitops','scrapeR','igraph','circlize','enrichR','readr','dplyr','gplots'))"


# RUN git clone https://github.com/sig2lead/sig2lead_beta.git

# RUN cd /sig2lead_beta/Sig2Lead_beta


#copy the app to the image
# COPY Sig2Lead_Backend.Rproj /srv/shiny-server/
# COPY ui.R /srv/shiny-server/
# COPY sample2.txt /srv/shiny-server/
# COPY LINCSCompounds.csv /srv/shiny-server/
# COPY lincs_fps.RData /srv/shiny-server/
# COPY minSim_apfp_RObjects.Rdata /srv/shiny-server/
# COPY ui.R /srv/shiny-server/
# COPY server.R /srv/shiny-server/
# COPY lib /srv/shiny-server/lib
# COPY similarity_search /srv/shiny-server/similarity_search

#select port
# EXPOSE 3838

#allow permission
# RUN sudo chown -R shiny:shiny /srv/shiny-server /var/lib/shiny-server /usr/bin

# RUN dos2unix /usr/bin/shiny-server.sh && apt-get --purge remove -y dos2unix && rm -rf /var/lib/apt/lists/*

#run app
# CMD ["/usr/bin/shiny-server.sh"]
