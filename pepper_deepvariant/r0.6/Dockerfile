FROM kishwars/pepper_deepvariant:r0.6

# for easy upgrade later. ARG variables only persist during image build
ARG PEPPERVER=0.6

# bugfix
RUN sed -i '238s/bam/fasta/' /usr/local/lib/python3.8/dist-packages/run_pepper_margin_deepvariant/run_pepper_margin_deepvariant.py && \
    sed -i '140s/^/# /' /usr/local/lib/python3.8/dist-packages/run_pepper_margin_deepvariant/run_pepper_margin_deepvariant.py