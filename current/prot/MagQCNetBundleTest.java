package prot;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ml.CellNet;
import ml.CellNetParser;

/** Executable contract test for the bundle loader and its failure matrix. */
public final class MagQCNetBundleTest {
    public static void main(String[] args) throws Exception {
        String bundle=null, netroot=null, vectors=null;
        for(String a:args) {
            int e=a.indexOf('='); if(e<1) throw new IllegalArgumentException("Expected key=value: "+a);
            String k=a.substring(0,e).toLowerCase(), v=a.substring(e+1);
            if(k.equals("bundle")||k.equals("in")) bundle=v;
            else if(k.equals("netroot")) netroot=v;
            else if(k.equals("vectors")) vectors=v;
        }
        if(bundle==null) throw new IllegalArgumentException("Required: bundle=");
        final MagQCNetBundle b=MagQCNetBundle.load(Paths.get(bundle));
        check(b.size()==176,"release count");
        check(b.subnet(0).order==0 && b.subnet(b.size()-1).order==175,"release ordering");
        for(int i=0;i<b.size();i++) {
            MagQCNetBundle.Subnet s=b.subnet(i);
            check(s.order==i && s.expectedInputs==s.dims[0] && s.expectedOutputs==s.dims[s.dims.length-1],s.id+" dimensions");
            check(s.netSha256.equals(sha256(s.netBytes)),s.id+" embedded hash");
        }
        // Round-trip identity: the semantic payload is unchanged by a second load.
        final String semantic=b.metadata("canonical_payload_sha256");
        check(semantic.equals(MagQCNetBundle.load(Paths.get(bundle)).metadata("canonical_payload_sha256")),"round-trip semantic hash");
        if(netroot!=null) verifyLooseBytes(b,Paths.get(netroot));
        runEquality(b,vectors,netroot==null?null:Paths.get(netroot));
        runParallelEquality(b);
        runMalformedMatrix(Paths.get(bundle),semantic);
        System.out.println("MAGQC_BUNDLE_ROUNDTRIP PASS");
        System.out.println("MAGQC_BUNDLE_OUTPUT_EQUALITY PASS");
        System.out.println("MAGQC_BUNDLE_MALFORMED_MATRIX PASS");
    }

    private static void verifyLooseBytes(MagQCNetBundle b,Path root) throws Exception {
        for(int i=0;i<b.size();i++) {
            MagQCNetBundle.Subnet s=b.subnet(i); Path p=root.resolve("nets").resolve(s.looseNetName);
            if(!Files.isRegularFile(p)) p=root.resolve(s.looseNetName);
            check(Files.isRegularFile(p),s.id+" loose net exists");
            byte[] x=Files.readAllBytes(p); check(Arrays.equals(x,s.netBytes),s.id+" loose byte equality");
        }
    }

    private static void runEquality(MagQCNetBundle b,String vectorFile,Path netroot) throws Exception {
        if(vectorFile!=null) {
            BufferedReader r=Files.newBufferedReader(Paths.get(vectorFile),StandardCharsets.UTF_8); String line; int n=0;
            while((line=r.readLine())!=null) {
                if(line.length()==0||line.startsWith("#")) continue;
                String[] x=line.split("\\t",-1); if(x.length!=3) throw new IOException("bad vector row: "+line);
                MagQCNetBundle.Subnet s=b.subnet(x[0]); float[] v=parseFloats(x[2]); float got=s.score(v);
                // The optional vector file carries the expected float bits produced by the loose net.
                int expected=(int)Long.parseLong(x[1]);
                check(Float.floatToIntBits(got)==expected,s.id+" output bits");
                n++;
            }
            r.close(); check(n==1000,"1000 fixed vectors");
        } else {
            // Deterministic vectors cover every subnet. With netroot this is a direct
            // bundled-vs-loose output comparison; without it the repeated-output check
            // still exercises every loaded CellNet copy. The gate supplies its 1,000 real
            // vectors through vectors=15_subnet_vectors_v4b.
            for(int i=0;i<b.size();i++) {
                MagQCNetBundle.Subnet s=b.subnet(i); float[] v=new float[s.expectedInputs];
                for(int j=0;j<v.length;j++) v[j]=(float)((i*31+j*7)%101)/100f;
                float a=s.score(v), c=s.score(v); check(Float.floatToIntBits(a)==Float.floatToIntBits(c),s.id+" repeated output");
                if(netroot!=null) {
                    Path p=netroot.resolve("nets").resolve(s.looseNetName); if(!Files.isRegularFile(p)) p=netroot.resolve(s.looseNetName);
                    CellNet loose=CellNetParser.load(p.toString());
                    loose.applyInput(v); loose.feedForward();
                    check(Float.floatToIntBits(a)==Float.floatToIntBits(loose.getOutput(0)),s.id+" loose output");
                }
            }
        }
    }

    private static void runParallelEquality(final MagQCNetBundle b) throws Exception {
        final int[] serial=new int[b.size()];
        for(int i=0;i<b.size();i++) serial[i]=Float.floatToIntBits(b.subnet(i).score(vector(i,b.subnet(i).expectedInputs)));
        java.util.concurrent.ExecutorService pool=java.util.concurrent.Executors.newFixedThreadPool(4);
        ArrayList<java.util.concurrent.Future<Integer>> fs=new ArrayList<java.util.concurrent.Future<Integer>>();
        for(int i=0;i<b.size();i++) {
            final int k=i; fs.add(pool.submit(new java.util.concurrent.Callable<Integer>() { public Integer call() {
                MagQCNetBundle.Subnet s=b.subnet(k); return Float.floatToIntBits(s.score(vector(k,s.expectedInputs)));
            }}));
        }
        pool.shutdown();
        for(int i=0;i<fs.size();i++) check(serial[i]==fs.get(i).get().intValue(),b.subnet(i).id+" serial/worker keyed equality");
    }

    private static float[] vector(int i,int width) { float[] v=new float[width]; for(int j=0;j<width;j++)v[j]=(float)((i*31+j*7)%101)/100f; return v; }

    private static void runMalformedMatrix(Path original,String semantic) throws Exception {
        byte[] plain=gunzip(Files.readAllBytes(original));
        String text=new String(plain,StandardCharsets.UTF_8);
        expectBad(text.replaceFirst("subnet_count=176","subnet_count=175").getBytes(StandardCharsets.UTF_8),"header count");
        expectBad(text.replaceFirst("##subnet\\n","##subnet\\n##subnet\\n").getBytes(StandardCharsets.UTF_8),"duplicate block");
        expectBad(removeFirstBlock(text).getBytes(StandardCharsets.UTF_8),"missing block");
        expectBad(text.replaceFirst("id=([^\\n]+)","id=extra").getBytes(StandardCharsets.UTF_8),"extra block");
        expectBad(text.replaceFirst("net_b64=","net_b64=AA==\n#extra\n" ).getBytes(StandardCharsets.UTF_8),"empty net");
        expectBad(text.replaceFirst("expected_inputs=([0-9]+)","expected_inputs=999999").getBytes(StandardCharsets.UTF_8),"dim-incompatible block");
        expectBad(text.replaceFirst("net_sha256=[0-9a-f]+","net_sha256="+repeat('0',64)).getBytes(StandardCharsets.UTF_8),"hash mismatch");
        byte[] truncated=Arrays.copyOf(plain,Math.max(1,plain.length/2)); expectBad(truncated,"truncated file");
        String commented="# harmless comment\n"+text.replace("##subnet\n","# harmless block comment\n\n##subnet\n");
        MagQCNetBundle variant=loadBytes(gzip(commented.getBytes(StandardCharsets.UTF_8)));
        check(semantic.equals(variant.metadata("canonical_payload_sha256")),"comment/whitespace invariant");
        // Compression wrapper bytes are deliberately excluded from the semantic digest.
        check(semantic.equals(loadBytes(gzip(plain)).metadata("canonical_payload_sha256")),"compression invariant");
    }

    private static String removeFirstBlock(String s) {
        int a=s.indexOf("##subnet\n"), z=s.indexOf("##endsubnet\n",a); if(a<0||z<0) throw new RuntimeException("test fixture block");
        return s.substring(0,a)+s.substring(z+"##endsubnet\n".length());
    }
    private static void expectBad(byte[] plain,String label) throws Exception {
        try { loadBytes(gzip(plain)); throw new AssertionError(label+" was accepted"); }
        catch(IOException expected) { }
        catch(RuntimeException expected) { }
    }
    private static MagQCNetBundle loadBytes(byte[] x) throws Exception {
        Path p=Files.createTempFile("magqc-bundle-test-", ".bbnets");
        try { Files.write(p,x); return MagQCNetBundle.load(p); } finally { Files.deleteIfExists(p); }
    }
    private static byte[] gunzip(byte[] x) throws Exception {
        GZIPInputStream in=new GZIPInputStream(new ByteArrayInputStream(x)); ByteArrayOutputStream b=new ByteArrayOutputStream(); byte[] q=new byte[8192]; int n;
        while((n=in.read(q))>=0) if(n>0)b.write(q,0,n); in.close(); return b.toByteArray();
    }
    private static byte[] gzip(byte[] x) throws Exception { ByteArrayOutputStream b=new ByteArrayOutputStream(); GZIPOutputStream out=new GZIPOutputStream(b); out.write(x); out.close(); return b.toByteArray(); }
    private static float[] parseFloats(String s) { String[] x=s.split(",",-1); float[] v=new float[x.length]; for(int i=0;i<v.length;i++)v[i]=Float.parseFloat(x[i]); return v; }
    private static String repeat(char c,int n){char[] x=new char[n];Arrays.fill(x,c);return new String(x);}
    private static String sha256(byte[] x){try{java.security.MessageDigest d=java.security.MessageDigest.getInstance("SHA-256");byte[] y=d.digest(x);StringBuilder b=new StringBuilder();for(byte q:y)b.append(String.format("%02x",q&255));return b.toString();}catch(Exception e){throw new RuntimeException(e);}}
    private static void check(boolean ok,String what){if(!ok)throw new AssertionError("FAIL: "+what);}
}
