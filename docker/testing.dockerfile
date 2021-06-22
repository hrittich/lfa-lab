FROM lfa-lab-base

COPY container-main.sh /usr/local/bin
RUN chmod a+x /usr/local/bin/container-main.sh

RUN useradd -m tester

WORKDIR /home/tester
ENTRYPOINT ["/usr/local/bin/container-main.sh"]
