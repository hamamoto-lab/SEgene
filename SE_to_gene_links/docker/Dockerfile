FROM condaforge/miniforge3:24.9.2-0

WORKDIR /home/jovyan/work

COPY environment.yml pyproject.toml setup.py /home/jovyan/code/

COPY src/ /home/jovyan/code/src/

RUN mkdir -p /root/.ipython/profile_default/startup/
COPY docker/startup/00-imports.py /root/.ipython/profile_default/startup/

RUN conda install -n base -c conda-forge conda=24.11.0 && \
    conda env create -f /home/jovyan/code/environment.yml && \
    conda clean -afy

ENV PATH=/opt/conda/envs/se_gene/bin:$PATH
ENV PYTHONPATH=/home/jovyan/code

RUN pip install -e /home/jovyan/code

EXPOSE 8888

RUN mkdir -p /home/jovyan/data


CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--allow-root", "--NotebookApp.token=''"]