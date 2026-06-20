package stream.bam;

import java.io.FileInputStream;

/**
 * Minimal test to exercise closing a BGZF MT stream while
 * producer/worker threads may be blocked on full queues.
 */
public class TestEarlyClose {
    //Test harness: opens BgzfInputStreamMT, sleeps 150 ms to let producer/workers fill queues,
    //then calls close() and prints elapsed time. No assertion — hangs indefinitely if close() deadlocks
    //(the known MT early-close deadlock is exactly what this targets). main()-only.
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

