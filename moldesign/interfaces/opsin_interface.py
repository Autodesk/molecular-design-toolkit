import moldesign as mdt
import moldesign.logs
from moldesign import utils

def name_to_smiles(name,
                   image='opsin',
                   engine=None):

    command = 'opsin -osmi input.txt output.txt'

    # TODO: this boilerplate has to go
    engine = utils.if_not_none(engine, mdt.compute.default_engine)
    imagename = mdt.compute.get_image_path(image, engine)
    job = engine.launch(imagename,
                          command,
                          inputs={'input.txt': name + '\n'},
                          name="opsin, %s" % name)
    moldesign.logs.display(job, "opsin, %s" % name)
    job.wait()
    return job.get_output('output.txt').read().strip()
