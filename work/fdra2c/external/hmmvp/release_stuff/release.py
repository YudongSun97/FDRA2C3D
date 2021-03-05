__all__ = ['make_release']

import os, sys, glob, string, re, shutil

hm_version = '1.3'

hm_home = '/home/ambrad/code/hmmvp/'
ut_home = '/home/ambrad/code/util/'
base_dir = '/scratch/ambrad/release/'
rm_certificate = base_dir + 'rm.txt'
hm_hdr = """\
/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version {0}
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/
""".format(hm_version)

def make_release():
    rdir = base_dir + 'hmmvp_v' + hm_version + '/'
    os.system('touch ' + rm_certificate)
    rmrfstar_with_protection(rdir)

    dirs = ['util/include', 'hmmvp/include', 'src', 'matlab', 'bin', 'lib',
            'examples', 'tmp'];
    [os.makedirs(rdir + d) for d in dirs]

    copy_dir(ut_home + 'include', rdir + 'util/include', '*.?pp', hm_hdr)
    copy_dir(hm_home + 'include', rdir + 'hmmvp/include', '*.hpp', hm_hdr)
    copy_dir(hm_home + 'include', rdir + 'hmmvp/include', '*.h', hm_hdr)

    copy_dir(ut_home + 'src', rdir + 'src', '*.?pp', hm_hdr)
    copy_dir(hm_home + 'src', rdir + 'src', '*.?pp', hm_hdr)

    copy_dir(ut_home + 'matlab', rdir + 'matlab', 'MexUtil.?pp', hm_hdr)
    copy_dir(hm_home + 'matlab', rdir + 'matlab', '*.?pp', hm_hdr)
    copy_file(ut_home + 'matlab/kvf.m', rdir + 'matlab')
    copy_file(hm_home + 'matlab/hmmvp.m', rdir + 'matlab')

    copy_file(hm_home + 'release_stuff/ex.m', rdir + 'examples')
    copy_file(hm_home + 'release_stuff/mvp_omp.cpp', rdir + 'examples')
    copy_file(hm_home + 'release_stuff/mvp_mpi.cpp', rdir + 'examples')
    copy_file(hm_home + 'release_stuff/cmvp_omp.c', rdir + 'examples')
    copy_file(hm_home + 'release_stuff/fmvp_omp.f90', rdir + 'examples')
    
    copy_file(hm_home + 'release_stuff/make.m', rdir + 'matlab')
    copy_file(hm_home + 'release_stuff/Makefile', rdir)
    prepend_header(hm_home + 'readme.txt', rdir, hm_hdr)

    os.system('cd ' + base_dir + '; ' +
              'rm -f ' + 'hmmvp_v' + hm_version + '.zip; ' +
              'zip -r ' + 'hmmvp_v' + hm_version + '.zip ' + 'hmmvp_v' +
              hm_version + '/*;')

def copy_dir(src_dir, dst_dir, glob_str, hdr):
    fns = glob.glob(src_dir + '/' + glob_str)
    [prepend_header(fn, dst_dir, hdr) for fn in fns]

def copy_file(src_fn, dst_dir):
    fn = os.path.split(src_fn)[1]
    shutil.copyfile(src_fn, dst_dir + '/' + fn)

def prepend_header(src_fn, dst_dir, hdr):
    fn = os.path.split(src_fn)[1]
    fr = open(src_fn, 'r')
    fw = open(dst_dir + '/' + fn, 'w')
    fw.write(hdr + '\n')
    for ln in fr:
        fw.write(ln)
    fr.close()
    fw.close()

def rmrfstar_with_protection(dr):
    """rm -rf $dr/*, but only if $dr/$rm_certificate exists, and only after
    additionally prompting."""
    if not os.path.exists(rm_certificate):
        print 'Not running rm -rf ' + dr + os.sep + '*'
        return
    cmd = 'rm -rf ' + dr + os.sep + '*'
    if ask_yes(cmd):
        os.system(cmd)

def ask_yes(prompt):
    while True:
        reply = raw_input(prompt + ' ' + '[[Yy]/[Nn]]: ')
        if reply in ['Y', 'y', 'N', 'n']: break
    return reply in ['Y', 'y']

if __name__ == '__main__':
    make_release()
