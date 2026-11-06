"""Portable source checks plus real shell execution of the app wrapper."""
from pathlib import Path
import json,os,re,shutil,subprocess,sys,tempfile,unittest

REPO=Path(__file__).resolve().parents[2]
BASE='ghcr.io/jaeyoonsung/dnmbsuite@sha256:2702a0615bfadd6e98b8ac67b38f45d5db6cbe4957aa90a9d371e9ab62d1f47c'

class ContainerRecipe(unittest.TestCase):
 def test_base_and_original_interface(self):
  text=(REPO/'Dockerfile.dnmbsuite').read_text()
  directives=[line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#')]
  self.assertEqual(next(x for x in directives if x.startswith('FROM ')),'FROM '+BASE)
  self.assertFalse(any(re.match(r'(ENTRYPOINT|CMD|WORKDIR|ENV)\s',x) for x in directives))
  self.assertIn('538a43d60df049a4959cb8a2ca8d704a21b6dbd10c439384aedbdb18bd02b5b7',text)
  self.assertIn('env -u R_HOME',text)
  self.assertNotIn('apt-get',text)
  self.assertNotIn('COPY . ',text)
  self.assertIn('--dependency-library /opt/biotools/lib/R/library',text)
 def test_resources_and_context(self):
  conf=json.loads((REPO/'docker/dnmbsuite/engine.json').read_text())
  self.assertEqual(conf['reference_cache'],'/resources/rebase')
  self.assertEqual(conf['hmm_database'],'/resources/rm-hmm/RM_HMMs.hmm')
  self.assertEqual(conf['r_libraries'][0],'/opt/defenseviz/native-r/library')
  self.assertEqual(conf['blastdbcmd'],'/opt/biotools/bin/blastdbcmd')
  ignore=(REPO/'Dockerfile.dnmbsuite.dockerignore').read_text().splitlines()
  self.assertEqual(ignore[0],'**')
  self.assertIn('!scripts/install_native_runtime.py',ignore)
  self.assertFalse(any(x.startswith('!') and any(y in x.lower() for y in ('unlazy','data/','fonts','output')) for x in ignore))
  self.assertTrue((REPO/'docker/dnmbsuite/CATH-LICENSE.txt').read_text().startswith('                    GNU GENERAL PUBLIC LICENSE'))
 def test_documented_commands(self):
  text=(REPO/'docs/dnmbsuite-container.md').read_text()
  for expected in ('Dockerfile.dnmbsuite','--platform linux/amd64','dst=/resources,readonly','dst=/input,readonly','--evaluation full','--entrypoint /opt/defenseviz/bin/defenseviz','--out /output/run001','read-only','Arial is not redistributed'):
   self.assertIn(expected,text)
 @unittest.skipUnless(shutil.which('sh'),'POSIX shell needed for container wrapper')
 def test_wrapper_preserves_environment_arguments_and_exit(self):
  with tempfile.TemporaryDirectory() as directory:
   tmp=Path(directory); fake=tmp/'cli.py'; wrapper=tmp/'wrapper'; fontdir=tmp/'mounted fonts';fontdir.mkdir()
   fake.write_text('import json,os,sys\nprint(json.dumps({"args":sys.argv[1:],"env":dict(os.environ)}))\nsys.exit(int(os.environ.get("TEST_EXIT","0")))\n')
   shim=tmp/'cli';shim.write_text('#!/bin/sh\nexec '+shlex_quote(sys.executable)+' '+shlex_quote(str(fake))+' "$@"\n');shim.chmod(0o755)
   fc=tmp/'fc-cache';fc.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$FONT_CALL"\n');fc.chmod(0o755)
   wrapper.write_text((REPO/'docker/dnmbsuite/defenseviz').read_text().replace('/opt/defenseviz/venv/bin/defenseviz',str(shim)).replace('/usr/local/share/fonts/defenseviz',shlex_quote(str(fontdir))))
   env=dict(os.environ,PATH=str(tmp)+os.pathsep+os.environ['PATH'],R_HOME='/wrong',R_LIBS='/unchanged',MPLCONFIGDIR=str(tmp/'mpl'),XDG_CACHE_HOME=str(tmp/'cache'),FONT_CALL=str(tmp/'font_call'))
   env.pop('DEFENSEVIZ_ENGINE_CONFIG',None)
   argv=['sh',str(wrapper),'from-suite','/data/path with spaces','--out','/output/a+b']
   result=subprocess.run(argv,capture_output=True,text=True,env=env);self.assertEqual(result.returncode,0,result.stderr)
   row=json.loads(result.stdout);self.assertEqual(row['args'],argv[2:]);self.assertNotIn('R_HOME',row['env'])
   for key in ('PATH','R_LIBS','MPLCONFIGDIR','XDG_CACHE_HOME'):self.assertEqual(row['env'][key],env[key])
   self.assertEqual(row['env']['DEFENSEVIZ_ENGINE_CONFIG'],'/opt/defenseviz/container/engine.json')
   self.assertEqual((tmp/'font_call').read_text().splitlines(),['-f',str(fontdir)])
   env.update(DEFENSEVIZ_ENGINE_CONFIG='/custom config.json',TEST_EXIT='17')
   result=subprocess.run(argv,capture_output=True,text=True,env=env);self.assertEqual(result.returncode,17)
   self.assertEqual(json.loads(result.stdout)['env']['DEFENSEVIZ_ENGINE_CONFIG'],'/custom config.json')

def shlex_quote(s):
 import shlex
 return shlex.quote(s)

if __name__=='__main__':unittest.main()
