FROM lfa-lab-base

RUN useradd -m tester
WORKDIR /home/tester

USER tester
#COPY --chown=tester:tester repository.git .
COPY --chown=tester:tester run-test-script.sh .
RUN chmod u+x run-test-script.sh

ENTRYPOINT ["/home/tester/run-test-script.sh"]

