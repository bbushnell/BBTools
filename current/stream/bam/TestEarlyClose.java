package stream.bam;

import java.io.FileInputStream;

/**
 * Minimal test to exercise closing a BGZF MT stream while
 * producer/worker threads may be blocked on full queues.
 */
public class TestEarlyClose {
    public static void main(String[] args) throws Exception {
        if(args.length<1){
            System.err.println("Usage: java stream.bam.TestEarlyClose <bgzf-file> [threads]");
            System.exit(1);
        }
        final String path=args[0];
        final int threads=(args.length>1? Integer.parseInt(args[1]) : 1);

        FileInputStream fis=new FileInputStream(path);
        BgzfInputStreamMT bgzf=new BgzfInputStreamMT(fis, threads);
        // Do not read; let producer/workers run a bit
        Thread.sleep(150);
        long t0=System.currentTimeMillis();
        bgzf.close();
        fis.close();
        long dt=System.currentTimeMillis()-t0;
        System.out.println("Closed cleanly in "+dt+" ms");
    }
}

