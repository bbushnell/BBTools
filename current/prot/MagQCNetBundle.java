package prot;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.BufferedOutputStream;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ml.CellNet;
import ml.CellNetParser;

/**
 * Self-contained, deterministic container for the MAG-QC subnet model set.
 *
 * <p>The outer file is always gzip-compressed UTF-8 text, regardless of its
 * suffix.  Net payloads are base64 encoded so CellNet text is preserved byte
 * for byte and cannot be confused with bundle records.  The loader constructs
 * and validates an unpublished object first; callers only receive it after all
 * header, manifest, digest, block, and CellNet gates pass.</n+ */
public final class MagQCNetBundle {

    private static final Charset UTF8 = StandardCharsets.UTF_8;
    private static final String MAGIC = "MAGQC_BBNETS_V1";
    private static final String HASH_ALGORITHM = "SHA-256";
    private static final int WRAP = 76;
    private static final String[] DOMAIN_VOCAB = {
        "bacteria", "archaea", "fungi", "plant", "animal", "protist", "virus", "other"
    };
    private static final String[] CONTEXT_ORDER = {
        "size", "gc", "coding", "genes", "l2rich", "glen", "glenstd", "hh", "caga",
        "domain_one_hot[8]"
    };
    /** Frozen four-output subnet contract (Part II.B), verbatim from the accepted
     * `two_stage_train.py` FOUR_OUTPUT_NAMES/FOUR_OUTPUT_UNITS constants, gate-proven at
     * bbtools-dev `264eaf8`. Order is fixed; a bundle entry's decoded names/units must equal
     * this exactly, not merely contain four non-empty strings. */
    private static final String[] FOUR_OUTPUT_NAMES = {
        "gene_completeness", "gene_contamination",
        "abs_residual_gene_completeness", "abs_residual_gene_contamination"
    };
    private static final String[] FOUR_OUTPUT_UNITS = {
        "fraction [0,1]",
        "fraction [0,1]",
        "absolute fraction-point error |pred-actual| for gene_completeness, same scale as column 0",
        "absolute fraction-point error |pred-actual| for gene_contamination, same scale as column 1"
    };
    private static final String SCHEMA2_OUTPUT_SEMANTICS =
        "per-subnet expected_outputs declares 1 (legacy) or 4 (learned); four-output entries carry "+
        "frozen output_names/output_units describing their own targets directly, with no derived "+
        "aggregator transform; legacy one-output entries remain consumed via the existing "+
        "[ratio,log_obs,log_pred,zero_flag] aggregator transform downstream, unchanged";
    private static final String FOUR_OUTPUT_CLAMP_RATIO =
        "no clamp applied by this API; output_units are target semantics";

    private final List<Subnet> subnets;
    private final Map<String, String> metadata;

    private MagQCNetBundle(List<Subnet> subnets, Map<String, String> metadata) {
        this.subnets = Collections.unmodifiableList(new ArrayList<Subnet>(subnets));
        this.metadata = Collections.unmodifiableMap(new LinkedHashMap<String, String>(metadata));
    }

    /** Inference-only CellNet subtype (root's prototype, `results/subnet_bundle_instance_dispatch_probe_20260909/`,
     * independently reviewed by Eru; released via `mag-qc/plans/THREAD_SAFE_BUNDLE_INFERENCE_20260909.md`
     * to replace the rejected static-DENSE-reset attempt, 2026-09-09). Constructs fresh,
     * structurally-empty cells via the public CellNet constructor, then copies weight/bias/input
     * state from the master net via setFrom(master,false) -- which does NOT call
     * makeWeightMatrices() (unlike copy()), so construction never touches the shared
     * ml.CellNet.DENSE global. feedForward() is overridden to dispatch directly on this
     * instance's OWN final denseMode flag -- this protects the DIRECT
     * netForCurrentThread().feedForward() call path MagQCVectorMaker already uses today, not
     * just the score()/scoreAll() wrappers, and involves no shared mutable state at inference
     * time, so no lock and no race between concurrently-scored different-format subnets.
     * Verified directly against BBTools source (2026-09-09): both the scalar
     * (ml.CellNet.feedForwardDense -> simd.Vector.feedForwardDense) and SIMD
     * (simd.SIMD.feedForward, reached when simdFF/Shared.SIMD/Shared.SIMD_FEED_FORWARD are all
     * true and width>=MINLEN32) dense code paths read only per-cell `weights`/`bias` -- exactly
     * what setFrom() populates -- never depending on makeWeightMatrices()'s separately-built
     * matrices; this construction is correct regardless of those flags' runtime values.
     * Inference only -- not for training, reserialization, or generic copy workflows. */
    private static final class InferenceNet extends CellNet {
        private final boolean denseMode;
        InferenceNet(CellNet master, int[] dims, long seed, float density, float density1,
                int blockSize, boolean denseMode) {
            super(dims, seed, density, density1, blockSize, new ArrayList<String>());
            this.denseMode=denseMode;
            setFrom(master, false);
        }
        @Override public float feedForward() { return denseMode ? feedForwardDense() : feedForwardSparse(); }
    }

    /** One immutable metadata record and a read-only master CellNet. */
    public static final class Subnet {
        public final String id, type, name, looseNetName, looseSubsetName;
        public final int order, numObs, expectedInputs, expectedOutputs, seed;
        public final int[] familyRanks, dims;
        public final byte[] netBytes;
        public final String netSha256, subsetSha256;
        /** Null for legacy (expectedOutputs==1) entries; frozen FOUR_OUTPUT_NAMES/UNITS
         * contract for expectedOutputs==4 entries (validated exactly against that contract
         * at load time, not merely checked for presence). Private and array-typed: exposed only
         * through outputNames()/outputUnits(), which return a fresh defensive copy each call, so
         * a caller cannot mutate this validated object's contract in place after construction
         * (a `final` array field is not itself immutable -- only the reference is). */
        private final String[] outputNames, outputUnits;
        /** This entry's OWN dense/sparse-ness, and its OWN density/density1/blockSize, all read
         * from its own raw net bytes at construction time -- NOT the shared global
         * ml.CellNet.DENSE flag (dense) and not hardcoded/prototype values (the other three).
         * Used only to construct a fresh InferenceNet per thread; never mutated after
         * construction. See parseDenseFlag()/parseDensities()/parseBlockSize(). */
        private final boolean dense;
        private final float density, density1;
        private final int blockSize;
        /** A SEPARATE clone of the dims array used only for InferenceNet construction, independent
         * of the public `dims` field. `dims` is `public final` -- the reference cannot be
         * reassigned, but the ARRAY IT POINTS TO is still mutable in place (an external caller
         * holding this Subnet could do `subnet.dims[0]=...`), and netForCurrentThread() calls
         * new InferenceNet(...) fresh on every cold ThreadLocal miss, so a later mutation of the
         * public field would silently feed a corrupted shape into a FUTURE construction (even
         * though CellNet's own constructor also clones its dims_ parameter -- ml/CellNet.java:75 --
         * that only protects an ALREADY-CONSTRUCTED instance, not the next one built from this
         * Subnet). Root's review, 2026-09-09: verified real by reading CellNet's constructor
         * directly, not assumed. */
        private final int[] executionDims;
        private final CellNet master;
        private final ThreadLocal<CellNet> local = new ThreadLocal<CellNet>();

        private Subnet(String id, String type, String name, int order, int numObs,
                int expectedInputs, int expectedOutputs, int seed, int[] familyRanks,
                int[] dims, byte[] netBytes, String netSha256, String subsetSha256,
                String looseNetName, String looseSubsetName, CellNet master,
                String[] outputNames, String[] outputUnits, boolean dense,
                float density, float density1, int blockSize) {
            this.id=id; this.type=type; this.name=name; this.order=order; this.numObs=numObs;
            this.expectedInputs=expectedInputs; this.expectedOutputs=expectedOutputs; this.seed=seed;
            this.familyRanks=familyRanks==null ? null : familyRanks.clone(); this.dims=dims.clone();
            this.executionDims=dims.clone(); // independent clone -- see field javadoc
            this.netBytes=netBytes.clone(); this.netSha256=netSha256; this.subsetSha256=subsetSha256;
            this.looseNetName=looseNetName; this.looseSubsetName=looseSubsetName; this.master=master;
            this.outputNames=outputNames==null ? null : outputNames.clone();
            this.outputUnits=outputUnits==null ? null : outputUnits.clone();
            this.dense=dense; this.density=density; this.density1=density1; this.blockSize=blockSize;
        }

        /** Fresh defensive copy each call; null for legacy (expectedOutputs==1) entries. */
        public String[] outputNames() { return outputNames==null ? null : outputNames.clone(); }
        /** Fresh defensive copy each call; null for legacy (expectedOutputs==1) entries. */
        public String[] outputUnits() { return outputUnits==null ? null : outputUnits.clone(); }

        /** Returns this worker's lazy private copy; no mutable master is exposed. A fresh
         * InferenceNet (not master.copy(false)) so construction never touches the shared
         * ml.CellNet.DENSE global, and this entry's own dense/sparse dispatch is baked into the
         * returned object's feedForward() override -- correct for direct virtual feedForward()
         * calls (as MagQCVectorMaker already makes) as well as the score()/scoreAll() wrappers. */
        public CellNet netForCurrentThread() {
            CellNet n=local.get();
            if(n==null) { n=new InferenceNet(master,executionDims,seed,density,density1,blockSize,dense); local.set(n); }
            return n;
        }

        /** Applies one vector to the current worker's copy and returns output zero. Legacy,
         * single-output entries only -- throws UNCONDITIONALLY (not via assert, so this check
         * fires even without -ea) for any other output count, so a caller cannot silently read
         * output zero of a multi-output net and mistake it for the whole prediction. Use
         * scoreAll() for expectedOutputs>1. */
        public float score(float[] input) {
            if(expectedOutputs!=1) {
                throw new IllegalStateException(id+": score() is single-output only (this entry "+
                        "declares "+expectedOutputs+" outputs) -- use scoreAll() to avoid silently "+
                        "reading output zero of a multi-output prediction");
            }
            if(input==null || input.length!=expectedInputs) {
                throw new IllegalArgumentException(id+": expected "+expectedInputs+" inputs, got "+
                        (input==null ? "null" : input.length));
            }
            final CellNet n=netForCurrentThread();
            n.applyInput(input); n.feedForward();
            return n.getOutput(0);
        }

        /** Applies one vector to the current worker's copy and returns every output, in order,
         * as an independent snapshot array (CellNet.getOutput() already allocates and copies
         * per cell -- confirmed by reading ml.CellNet.java lines 1035-1043 -- so no further
         * defensive copy is needed here). Works for any output count, including legacy
         * single-output entries (returns a length-1 array). */
        public float[] scoreAll(float[] input) {
            if(input==null || input.length!=expectedInputs) {
                throw new IllegalArgumentException(id+": expected "+expectedInputs+" inputs, got "+
                        (input==null ? "null" : input.length));
            }
            final CellNet n=netForCurrentThread();
            n.applyInput(input); n.feedForward();
            return n.getOutput();
        }
    }

    public int size() { return subnets.size(); }
    public Subnet subnet(int index) { return subnets.get(index); }
    public Subnet subnet(String id) {
        for(Subnet s : subnets) if(s.id.equals(id)) return s;
        throw new IllegalArgumentException("Unknown subnet: "+id);
    }
    public String metadata(String key) { return metadata.get(key); }

    /** Fails closed unless the consumer's family-list bytes are exactly the pinned release bytes. */
    public void requireFamilyList(java.nio.file.Path familyList) throws IOException {
        final byte[] bytes=Files.readAllBytes(familyList);
        if(!sha256(bytes).equals(metadata("family_list_sha256"))) {
            throw new IOException("bundle family_list_sha256 does not match "+familyList);
        }
    }

    public static void main(String[] args) throws Exception {
        if(args.length<1) throw new IllegalArgumentException("Usage: pack|verify|list key=value ...");
        final String command=args[0].toLowerCase();
        final Map<String,String> a=parseArgs(args, 1);
        if("pack".equals(command)) pack(a);
        else if("verify".equals(command)) verify(a);
        else if("list".equals(command)) list(a);
        else throw new IllegalArgumentException("Unknown command: "+args[0]);
    }

    /** Packs the selected best net for each release-manifest subnet and atomically publishes it. */
    public static void pack(Map<String,String> a) throws Exception {
        final Path agg=requiredPath(a,"aggmanifest"), release=requiredPath(a,"subnetmanifest");
        final Path familyList=requiredPath(a,"familylist"), taxpgm=requiredPath(a,"taxpgm");
        final Path out=requiredPath(a,"out");
        final Path root=requiredPath(a,"netroot");
        final Path subsetRoot=pathOr(a,"subsetroot",root);
        final byte[] releaseBytes=Files.readAllBytes(release);
        final List<ReleaseRow> rows=readRelease(release);
        final Map<String,AggRow> aggRows=readAgg(agg);
        if(rows.size()!=176) throw new IOException("release manifest must contain exactly 176 rows, got "+rows.size());
        if(aggRows.size()!=rows.size()) throw new IOException("agg/release row count mismatch");
        final ArrayList<Subnet> records=new ArrayList<Subnet>(rows.size());
        final TreeSet<Integer> familyUniverse=new TreeSet<Integer>();
        for(int i=0; i<rows.size(); i++) {
            final ReleaseRow rr=rows.get(i);
            final AggRow ar=aggRows.get(rr.id);
            if(ar==null) throw new IOException("release subnet missing from aggmanifest: "+rr.id);
            if(!rr.type.equals(ar.type)) throw new IOException(rr.id+": type mismatch release="+rr.type+" agg="+ar.type);
            final Path net=resolve(root,ar.netPath);
            final byte[] netBytes=Files.readAllBytes(net);
            final CellNet cell=parseNet(netBytes, rr.id);
            final int[] dims=parseDims(netBytes,rr.id);
            final int expectedInputs=ar.expectedInputs>=0 ? ar.expectedInputs : cell.numInputs();
            if(rr.expectedOutputs!=1 && rr.expectedOutputs!=4) {
                throw new IOException(rr.id+": unsupported declared output count "+rr.expectedOutputs);
            }
            if(cell.numOutputs()!=rr.expectedOutputs) throw new IOException(rr.id+": manifest declares "+
                    rr.expectedOutputs+" output(s), net has "+cell.numOutputs());
            if(dims[dims.length-1]!=rr.expectedOutputs) throw new IOException(rr.id+": dims last entry "+
                    dims[dims.length-1]+" != declared output count "+rr.expectedOutputs);
            if(cell.numInputs()!=expectedInputs) throw new IOException(rr.id+": input width mismatch agg="+
                    expectedInputs+" net="+cell.numInputs());
            final int[] ranks;
            final Path subset;
            final String subsetName;
            if("ncrna".equals(rr.type)) {
                if(!"-".equals(ar.subsetPath) || !"-".equals(rr.subsetPath)) throw new IOException(rr.id+": ncrna subset must be -");
                ranks=null; subset=null; subsetName="-";
                if(ar.numObs!=5) throw new IOException(rr.id+": ncrna observation count must be 5");
            } else {
                subset=resolve(subsetRoot,rr.subsetPath);
                ranks=readRanks(subset,rr.id);
                subsetName=subset.getFileName().toString();
                if(ranks.length!=ar.numObs) throw new IOException(rr.id+": subset rank count="+ranks.length+" agg numObs="+ar.numObs);
                for(int r:ranks) { if(r<0) throw new IOException(rr.id+": negative family rank "+r); familyUniverse.add(r); }
            }
            final int seed=parseSeed(netBytes,rr.id);
            final boolean dense=parseDenseFlag(netBytes,rr.id);
            final float[] densities=parseDensities(netBytes,rr.id);
            final int blockSize=parseBlockSize(netBytes,rr.id);
            final String[] outNames=rr.expectedOutputs==4 ? FOUR_OUTPUT_NAMES : null;
            final String[] outUnits=rr.expectedOutputs==4 ? FOUR_OUTPUT_UNITS : null;
            records.add(new Subnet(rr.id,rr.type,rr.id,i,ar.numObs,expectedInputs,rr.expectedOutputs,seed,
                    ranks,dims,netBytes,sha256(netBytes),subset==null?"-":sha256(Files.readAllBytes(subset)),
                    net.getFileName().toString(),subsetName,cell,outNames,outUnits,dense,
                    densities[0],densities[1],blockSize));
        }
        final LinkedHashMap<String,String> md=baseMetadata(a,releaseBytes,rows,familyUniverse,familyList,taxpgm);
        md.put("canonical_payload_sha256",canonicalHash(md,records));
        writeAtomic(out,md,records);
        // Re-open the unpublished temporary-equivalent output. Always via the multi-output-
        // capable loader: pack() must be able to verify its own schema2 output too, and this
        // loader also accepts schema1 (pure legacy), so one call covers both cases.
        final MagQCNetBundle checked=loadMultiOutput(out);
        if(!md.get("canonical_payload_sha256").equals(checked.metadata("canonical_payload_sha256"))) {
            throw new IOException("post-publish canonical hash changed");
        }
        System.out.println("MAGQC_BBNETS_PACK PASS count="+checked.size()+" out="+out+
                " semantic_sha256="+checked.metadata("canonical_payload_sha256"));
    }

    private static LinkedHashMap<String,String> baseMetadata(Map<String,String> a, byte[] releaseBytes,
            List<ReleaseRow> rows, TreeSet<Integer> familyUniverse, Path familyList, Path taxpgm) throws Exception {
        boolean anyFourOutput=false;
        for(ReleaseRow r:rows) { if(r.expectedOutputs==4) { anyFourOutput=true; break; } }
        final LinkedHashMap<String,String> md=new LinkedHashMap<String,String>();
        md.put("schema_version",anyFourOutput?"2":"1");
        md.put("hash_algorithm",HASH_ALGORITHM);
        md.put("subnet_count",Integer.toString(rows.size()));
        md.put("release_manifest_sha256",sha256(releaseBytes));
        md.put("release_manifest_values_b64",b64(releaseBytes));
        final byte[] familyBytes=Files.readAllBytes(familyList);
        final String familyValues=readFamilyRepIds(familyBytes);
        md.put("family_list_definition","familylist_v4b.tsv ordered rep_id values; rank is the zero-based row index");
        md.put("family_list_values",familyValues);
        md.put("family_list_sha256",sha256(familyBytes));
        md.put("family_list_rank_universe",join(familyUniverse.stream().mapToInt(Integer::intValue).toArray(),","));
        md.put("vector_layout","[ordered subnet observations] + [phylum one-hot] + [9 shared context scalars] + [8 domain one-hot]");
        md.put("taxonomy_layout","phylum one-hot is release-consumer ordered vocabulary; domain one-hot uses frozen source mapping");
        final byte[] taxBytes=Files.readAllBytes(taxpgm);
        md.put("taxonomy_source_sha256",sha256(taxBytes));
        md.put("phylum_vocab_values",readPhylumVocabulary(taxBytes));
        md.put("domain_vocab_values",join(DOMAIN_VOCAB,","));
        md.put("context_order_values",join(CONTEXT_ORDER,","));
        md.put("context_constants","NOT_FROZEN_v4b");
        md.put("context_normalization_values","NOT_FROZEN_v4b; source computes corpus scales live after B2; no values are embedded or recomputed here");
        md.put("presence_copy_transforms","ratio=N/(1+N);raw=min(N,32)/32;log=log2(1+N)/log2(65);two=presence plus min(N-1,16)/16;norm=N/avgCopyWhenPresent");
        md.put("output_semantics",anyFourOutput ? SCHEMA2_OUTPUT_SEMANTICS :
                "one regression output; aggregator observations are [ratio,log_obs,log_pred,zero_flag]; ratio clamp=2.0; ncRNA observations=5");
        md.put("input_root_recorded",pathOr(a,"netroot",Paths.get(".")).toAbsolutePath().normalize().toString());
        md.put("subset_root_recorded",pathOr(a,"subsetroot",pathOr(a,"netroot",Paths.get(".")).toAbsolutePath()).toAbsolutePath().normalize().toString());
        md.put("order_rule","release manifest order; input file order and worker count do not affect bytes");
        return md;
    }

    private static String readFamilyRepIds(byte[] bytes) throws Exception {
        StringBuilder b=new StringBuilder(); String[] lines=new String(bytes,UTF8).split("\\r?\\n",-1); int rank=0;
        for(String line:lines) { if(line.length()==0 || line.startsWith("#")) continue; String[] x=line.split("\\t",-1); if(x.length<2) throw new IOException("bad family list row: "+line);
            if(Integer.parseInt(x[0])!=rank++) throw new IOException("family list rank discontinuity at "+line); if(b.length()>0)b.append(','); b.append(x[1]); }
        if(rank==0) throw new IOException("empty family list"); return b.toString();
    }
    private static String readPhylumVocabulary(byte[] bytes) throws Exception {
        TreeSet<String> x=new TreeSet<String>(); String[] lines=new String(bytes,UTF8).split("\\r?\\n",-1);
        for(String line:lines) { if(line.length()==0 || line.startsWith("#")) continue; String[] f=line.split("\\t",-1); if(f.length<2) throw new IOException("bad taxpgm row: "+line); x.add(f[1]); }
        if(x.isEmpty()) throw new IOException("empty taxonomy vocabulary"); StringBuilder b=new StringBuilder(); for(String s:x){if(b.length()>0)b.append(',');b.append(s);} b.append(",other"); return b.toString();
    }

    private static void writeAtomic(Path out, Map<String,String> md, List<Subnet> records) throws Exception {
        final Path parent=out.toAbsolutePath().normalize().getParent();
        if(parent!=null) Files.createDirectories(parent);
        final Path tmp=Files.createTempFile(parent,".bbnets-", ".tmp");
        boolean published=false;
        try {
            writeBundle(tmp,md,records);
            loadMultiOutput(tmp); // all gates pass before publication; accepts schema1 or schema2
            try { Files.move(tmp,out,StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING); }
            catch(AtomicMoveNotSupportedException e) { Files.move(tmp,out,StandardCopyOption.REPLACE_EXISTING); }
            published=true;
        } finally { if(!published) Files.deleteIfExists(tmp); }
    }

    private static void writeBundle(Path out, Map<String,String> md, List<Subnet> records) throws Exception {
        final OutputStream fos=new BufferedOutputStream(Files.newOutputStream(out));
        final GZIPOutputStream gz=new GZIPOutputStream(fos);
        final BufferedWriter w=new BufferedWriter(new OutputStreamWriter(gz,UTF8));
        w.write(MAGIC); w.write('\n');
        for(Map.Entry<String,String> e:md.entrySet()) {
            if("release_manifest_values_b64".equals(e.getKey())) writeWrapped(w,e.getKey(),e.getValue());
            else { w.write(e.getKey()); w.write('='); w.write(e.getValue()); w.write('\n'); }
        }
        for(Subnet s:records) {
            w.write("##subnet\n");
            field(w,"id",s.id); field(w,"name",s.name); field(w,"type",s.type); field(w,"order",s.order);
            field(w,"num_obs",s.numObs); field(w,"expected_inputs",s.expectedInputs); field(w,"expected_outputs",s.expectedOutputs);
            field(w,"seed",s.seed); field(w,"dims",join(s.dims,","));
            field(w,"family_ranks",ranksToken(s.familyRanks));
            field(w,"ncrna_obs_definition",s.familyRanks==null?"5 ordered ncRNA observations":"-");
            if(s.expectedOutputs==4) {
                field(w,"output_names",b64Lines(s.outputNames()));
                field(w,"output_units",b64Lines(s.outputUnits()));
                field(w,"output_clamp_ratio",FOUR_OUTPUT_CLAMP_RATIO);
            } else {
                field(w,"output_units","completeness_or_contamination_regression");
                field(w,"output_clamp_ratio","[0,2] ratio clamp; one scalar output");
            }
            field(w,"loose_net_name",s.looseNetName); field(w,"loose_subset_name",s.looseSubsetName);
            field(w,"subset_sha256",s.subsetSha256); field(w,"net_sha256",s.netSha256); field(w,"net_bytes",s.netBytes.length);
            writeWrapped(w,"net_b64",b64(s.netBytes));
            w.write("##endsubnet\n");
        }
        w.flush(); w.close();
    }

    private static void field(BufferedWriter w,String k,Object v) throws IOException { w.write(k); w.write('='); w.write(String.valueOf(v)); w.write('\n'); }
    private static void writeWrapped(BufferedWriter w,String key,String value) throws IOException {
        int p=0; boolean first=true;
        while(p<value.length()) { int e=Math.min(p+WRAP,value.length()); w.write(first?key+"=":"+"); w.write(value.substring(p,e)); w.write('\n'); first=false; p=e; }
        if(first) w.write(key+"=\n");
    }

    /** Loads only fully validated, schema_version=1 (legacy, single-output-only) bundles.
     * Rejects a schema_version=2 bundle outright -- this is the default entry point every
     * current production caller (MVM, MagQCTool, MagQCCompositeBundle) already uses, so none of
     * them can silently load a bundle containing a four-output subnet and mistake output zero
     * for the whole prediction. Use loadMultiOutput() to accept schema_version=2. */
    public static MagQCNetBundle load(Path in) throws Exception { return loadInternal(in,false); }

    /** Loads a fully validated bundle of either schema_version -- the only way to accept a
     * bundle containing a four-output subnet. An explicit opt-in a caller must deliberately
     * call; also accepts pure-legacy (schema_version=1) bundles, so one method covers both. */
    public static MagQCNetBundle loadMultiOutput(Path in) throws Exception { return loadInternal(in,true); }

    private static MagQCNetBundle loadInternal(Path in, boolean allowMultiOutput) throws Exception {
        final ArrayList<String> lines=new ArrayList<String>();
        final InputStream raw=new BufferedInputStream(Files.newInputStream(in));
        final GZIPInputStream gz;
        try { gz=new GZIPInputStream(raw); }
        catch(IOException e) { try { raw.close(); } catch(IOException ignored) {} throw new IOException("not a gzip .bbnets: "+in,e); }
        final BufferedReader r=new BufferedReader(new InputStreamReader(gz,UTF8));
        for(String line=r.readLine(); line!=null; line=r.readLine()) lines.add(line);
        r.close();
        int magicLine=0;
        while(magicLine<lines.size() && (lines.get(magicLine).length()==0 || lines.get(magicLine).startsWith("#"))) magicLine++;
        if(magicLine>=lines.size() || !MAGIC.equals(lines.get(magicLine))) throw new IOException("bad bundle magic");
        final LinkedHashMap<String,String> md=new LinkedHashMap<String,String>();
        final ArrayList<SubnetText> blocks=new ArrayList<SubnetText>();
        int p=magicLine+1;
        while(p<lines.size() && !"##subnet".equals(lines.get(p))) {
            final String line=lines.get(p++);
            if(line.length()==0 || line.startsWith("#")) continue;
            final int eq=line.indexOf('=');
            if(eq<=0 || line.startsWith("+")) throw new IOException("malformed header line: "+line);
            final String key=line.substring(0,eq), value=line.substring(eq+1);
            if(md.containsKey(key)) throw new IOException("duplicate header key: "+key);
            if("release_manifest_values_b64".equals(key)) {
                StringBuilder b=new StringBuilder(value); while(p<lines.size() && lines.get(p).startsWith("+")) b.append(lines.get(p++).substring(1));
                md.put(key,b.toString());
            } else md.put(key,value);
        }
        while(p<lines.size()) {
            while(p<lines.size() && (lines.get(p).length()==0 ||
                    (lines.get(p).startsWith("#") && !"##subnet".equals(lines.get(p))))) p++;
            if(p>=lines.size()) break;
            if(!"##subnet".equals(lines.get(p++))) throw new IOException("expected ##subnet at line "+p);
            final LinkedHashMap<String,String> f=new LinkedHashMap<String,String>();
            String netB64=null;
            while(p<lines.size() && !"##endsubnet".equals(lines.get(p))) {
                final String line=lines.get(p++);
                if(line.length()==0 || line.startsWith("#")) continue;
                final int eq=line.indexOf('=');
                if(eq<=0) throw new IOException("malformed subnet field: "+line);
                final String key=line.substring(0,eq), value=line.substring(eq+1);
                if("net_b64".equals(key)) {
                    if(netB64!=null) throw new IOException("duplicate subnet field: net_b64");
                    StringBuilder b=new StringBuilder(value); while(p<lines.size() && lines.get(p).startsWith("+")) b.append(lines.get(p++).substring(1));
                    netB64=b.toString();
                } else if(f.put(key,value)!=null) throw new IOException("duplicate subnet field: "+key);
            }
            if(p>=lines.size()) throw new IOException("truncated subnet block");
            p++; // end marker
            if(netB64==null) throw new IOException("empty subnet net payload");
            blocks.add(new SubnetText(f,netB64));
        }
        return validate(md,blocks,allowMultiOutput);
    }

    private static MagQCNetBundle validate(Map<String,String> md,List<SubnetText> blocks,boolean allowMultiOutput) throws Exception {
        require(md,"schema_version"); require(md,"hash_algorithm"); require(md,"subnet_count");
        require(md,"release_manifest_sha256"); require(md,"release_manifest_values_b64");
        require(md,"family_list_sha256"); require(md,"family_list_values"); require(md,"canonical_payload_sha256");
        final String schemaVersion=md.get("schema_version");
        if((!"1".equals(schemaVersion) && !"2".equals(schemaVersion)) || !HASH_ALGORITHM.equals(md.get("hash_algorithm"))) {
            throw new IOException("unsupported schema/hash");
        }
        if("2".equals(schemaVersion) && !allowMultiOutput) {
            throw new IOException("schema_version=2 bundle requires loadMultiOutput(), not load()");
        }
        // Metadata being hash-bound only proves internal self-consistency (the stored canonical
        // hash matches the stored content) -- it does NOT prove the content is the one true
        // frozen contract, since a forger can change a value AND recompute the hash. Require the
        // exact literal independently (root's review, 2026-09-09).
        if("2".equals(schemaVersion) && !SCHEMA2_OUTPUT_SEMANTICS.equals(md.get("output_semantics"))) {
            throw new IOException("schema_version=2 bundle's output_semantics does not match the frozen contract");
        }
        final int count=parseInt(md,"subnet_count");
        if(count!=blocks.size()) throw new IOException("header count="+count+" block count="+blocks.size());
        final byte[] release=decode(md.get("release_manifest_values_b64"),"release manifest");
        if(!sha256(release).equals(md.get("release_manifest_sha256"))) throw new IOException("release manifest hash mismatch");
        final List<ReleaseRow> expected=parseReleaseBytes(release);
        if(expected.size()!=count) throw new IOException("embedded release manifest count mismatch");
        final ArrayList<Subnet> result=new ArrayList<Subnet>(count); final Set<String> seen=new HashSet<String>();
        for(int i=0;i<blocks.size();i++) {
            final SubnetText st=blocks.get(i); final Map<String,String> f=st.fields;
            final String id=req(f,"id");
            if(!safeKey(id) || !seen.add(id)) throw new IOException("duplicate/unsafe subnet id: "+id);
            final ReleaseRow rr=expected.get(i);
            if(!id.equals(rr.id) || !req(f,"type").equals(rr.type) || parseInt(f,"order")!=i) throw new IOException("release order/id/type mismatch at "+i);
            final int numObs=parseInt(f,"num_obs"), expIn=parseInt(f,"expected_inputs"), expOut=parseInt(f,"expected_outputs");
            if(expOut!=1 && expOut!=4) throw new IOException(id+": unsupported declared output count "+expOut);
            if(rr.expectedOutputs!=expOut) throw new IOException(id+": release manifest declares "+
                    rr.expectedOutputs+" output(s) but subnet block declares "+expOut);
            final String[] outNames, outUnits;
            if(expOut==4) {
                final String[] names=decodeLines(req(f,"output_names"),id+" output_names");
                final String[] units=decodeLines(req(f,"output_units"),id+" output_units");
                if(!Arrays.equals(names,FOUR_OUTPUT_NAMES)) throw new IOException(id+": output_names does not match the frozen four-output contract");
                if(!Arrays.equals(units,FOUR_OUTPUT_UNITS)) throw new IOException(id+": output_units does not match the frozen four-output contract");
                if(!FOUR_OUTPUT_CLAMP_RATIO.equals(req(f,"output_clamp_ratio"))) throw new IOException(id+": output_clamp_ratio does not match the frozen four-output contract");
                outNames=names; outUnits=units;
            } else {
                if(f.containsKey("output_names")) throw new IOException(id+": legacy (1-output) entry must not declare output_names");
                outNames=null; outUnits=null;
            }
            final int seed=parseInt(f,"seed"); final int[] dims=parseInts(req(f,"dims"),",",id+" dims");
            final String rankToken=req(f,"family_ranks");
            final int[] ranks="-".equals(rankToken)?null:("empty".equals(rankToken)?new int[0]:parseInts(rankToken,",",id+" family ranks"));
            if("ncrna".equals(rr.type) != (ranks==null)) throw new IOException(id+": ncrna/famset payload mismatch");
            if(ranks==null && numObs!=5) throw new IOException(id+": bad ncRNA observation definition");
            if(ranks!=null && ranks.length!=numObs) throw new IOException(id+": rank count mismatch");
            final byte[] net=decode(st.netB64,id+" net");
            final String hash=sha256(net);
            if(!hash.equals(req(f,"net_sha256"))) throw new IOException(id+": net hash mismatch");
            if(parseInt(f,"net_bytes")!=net.length) throw new IOException(id+": net byte count mismatch");
            final CellNet cell=parseNet(net,id);
            if(cell.numInputs()!=expIn || cell.numOutputs()!=expOut || !Arrays.equals(dims,parseDims(net,id))) throw new IOException(id+": incompatible CellNet dimensions");
            final int actualSeed=parseSeed(net,id); if(actualSeed!=seed) throw new IOException(id+": seed mismatch");
            final boolean dense=parseDenseFlag(net,id);
            final float[] densities=parseDensities(net,id);
            final int blockSize=parseBlockSize(net,id);
            result.add(new Subnet(id,req(f,"type"),req(f,"name"),i,numObs,expIn,expOut,seed,ranks,dims,net,hash,
                    req(f,"subset_sha256"),req(f,"loose_net_name"),req(f,"loose_subset_name"),cell,outNames,outUnits,dense,
                    densities[0],densities[1],blockSize));
        }
        if(seen.size()!=expected.size()) throw new IOException("missing or extra release subnet");
        // schema_version must agree with what the entries actually contain, not just be a label.
        boolean anyFourOutput=false;
        for(Subnet s:result) { if(s.expectedOutputs==4) { anyFourOutput=true; break; } }
        if("1".equals(schemaVersion) && anyFourOutput) throw new IOException("schema_version=1 bundle contains a four-output subnet");
        if("2".equals(schemaVersion) && !anyFourOutput) throw new IOException("schema_version=2 bundle contains no four-output subnet");
        final String canonical=canonicalHash(md,result);
        if(!canonical.equals(md.get("canonical_payload_sha256"))) throw new IOException("semantic payload hash mismatch");
        return new MagQCNetBundle(result,md);
    }

    public static void verify(Map<String,String> a) throws Exception {
        final Path in=requiredPath(a,"in");
        final boolean multi=parseMultioutputFlag(a.get("multioutput"));
        final MagQCNetBundle b= multi ? loadMultiOutput(in) : load(in);
        final String releasePath=a.get("subnetmanifest");
        if(releasePath!=null) {
            final byte[] bytes=Files.readAllBytes(Paths.get(releasePath));
            if(!sha256(bytes).equals(b.metadata("release_manifest_sha256")) || !Arrays.equals(bytes,decode(b.metadata("release_manifest_values_b64"),"manifest")))
                throw new IOException("external release manifest mismatch");
        }
        final String root=a.get("netroot");
        if(root!=null) for(Subnet s:b.subnets) {
            final Path loose=resolve(Paths.get(root),s.looseNetName);
            final byte[] bytes=Files.readAllBytes(loose);
            if(!sha256(bytes).equals(s.netSha256) || !Arrays.equals(bytes,s.netBytes)) throw new IOException(s.id+": loose net mismatch");
        }
        System.out.println("MAGQC_BBNETS_VERIFY PASS count="+b.size()+" semantic_sha256="+b.metadata("canonical_payload_sha256"));
    }

    public static void list(Map<String,String> a) throws Exception {
        final boolean multi=parseMultioutputFlag(a.get("multioutput"));
        final MagQCNetBundle b= multi ? loadMultiOutput(requiredPath(a,"in")) : load(requiredPath(a,"in"));
        System.out.println("count="+b.size()+" semantic_sha256="+b.metadata("canonical_payload_sha256"));
        for(Subnet s:b.subnets) System.out.println(s.order+"\t"+s.id+"\t"+s.type+"\t"+s.expectedInputs+"\t"+s.netSha256);
    }

    private static String canonicalHash(Map<String,String> md,List<Subnet> records) {
        final boolean v2="2".equals(md.get("schema_version"));
        final StringBuilder b=new StringBuilder();
        // Pure-legacy bundles (schema_version=1) use the exact V1 domain literal and line
        // format, byte for byte -- this is the whole backward-compatibility guarantee, verified
        // by construction: a bundle with zero four-output entries never enters the v2 branches
        // below, so its canonical hash is computed identically to before this change.
        b.append(v2 ? "MAGQC_BBNETS_CANONICAL_V2\n" : "MAGQC_BBNETS_CANONICAL_V1\n");
        for(Map.Entry<String,String> e:md.entrySet()) {
            // These are provenance-only locations, not semantic bundle content.
            // Packing the same bytes from another staging directory must produce
            // the same semantic digest.
            if("canonical_payload_sha256".equals(e.getKey()) || "input_root_recorded".equals(e.getKey()) ||
                    "subset_root_recorded".equals(e.getKey())) continue;
            b.append(e.getKey()).append('=').append(e.getValue()).append('\n');
        }
        for(Subnet s:records) {
            b.append("subnet\t").append(s.order).append('\t').append(s.id).append('\t').append(s.name).append('\t').append(s.type).append('\t')
                .append(s.numObs).append('\t').append(s.expectedInputs).append('\t').append(s.expectedOutputs).append('\t').append(s.seed).append('\t')
                .append(join(s.dims,",")).append('\t').append(ranksToken(s.familyRanks)).append('\t')
                .append(s.looseNetName).append('\t').append(s.looseSubsetName).append('\t').append(s.subsetSha256).append('\t')
                .append(s.netSha256).append('\t').append(s.netBytes.length);
            // Legacy entries' lines are otherwise byte-identical to V1 -- only four-output
            // entries append their names/units/clamp-ratio, and only those entries ever have
            // non-null outputNames/outputUnits (a pure-legacy bundle never reaches this branch
            // at all). output_clamp_ratio is included here (root's review, 2026-09-09): it was
            // previously neither validated nor canonicalized for four-output entries, so a
            // tampered value survived a recomputed-checksum forgery undetected.
            if(s.expectedOutputs==4) {
                b.append('\t').append(b64Lines(s.outputNames())).append('\t').append(b64Lines(s.outputUnits()))
                 .append('\t').append(FOUR_OUTPUT_CLAMP_RATIO);
            }
            b.append('\n');
        }
        return sha256(b.toString().getBytes(UTF8));
    }

    private static List<ReleaseRow> readRelease(Path p) throws Exception { return parseReleaseBytes(Files.readAllBytes(p)); }
    /** A row is 3 tab-separated fields (id, type, subsetPath) -- outputs=1 implicit, byte-
     * identical to every existing manifest file, no migration needed -- OR 4 fields, where the
     * 4th MUST be the exact literal "outputs=4"; any other 4th-field content is rejected
     * outright (fail closed on unknown flags, never silently ignored or guessed at). */
    private static List<ReleaseRow> parseReleaseBytes(byte[] bytes) throws Exception {
        final ArrayList<ReleaseRow> out=new ArrayList<ReleaseRow>();
        final String text=new String(bytes,UTF8); final String[] ls=text.split("\\r?\\n",-1);
        final Set<String> seen=new HashSet<String>();
        for(String line:ls) { if(line.length()==0 || line.startsWith("#")) continue; String[] x=line.split("\\t",-1);
            if(x.length!=3 && x.length!=4) throw new IOException("bad release row: "+line);
            if(!safeKey(x[0]) || !("famset".equals(x[1])||"ncrna".equals(x[1])) || !seen.add(x[0])) throw new IOException("bad release row: "+line);
            int expectedOutputs=1;
            if(x.length==4) { if(!"outputs=4".equals(x[3])) throw new IOException("bad release row flag: "+line); expectedOutputs=4; }
            out.add(new ReleaseRow(x[0],x[1],x[2],expectedOutputs)); }
        return out;
    }
    private static Map<String,AggRow> readAgg(Path p) throws Exception {
        final LinkedHashMap<String,AggRow> out=new LinkedHashMap<String,AggRow>();
        final String[] ls=new String(Files.readAllBytes(p),UTF8).split("\\r?\\n",-1);
        for(String line:ls) { if(line.length()==0 || line.startsWith("#")) continue; String[] x=line.split("\\t",-1); if(x.length!=6) throw new IOException("bad agg row: "+line);
            if(!safeKey(x[0]) || out.containsKey(x[0])) throw new IOException("bad/duplicate agg id: "+x[0]);
            out.put(x[0],new AggRow(x[0],"ncrna".equals(x[0])?"ncrna":"famset",Integer.parseInt(x[1]),Integer.parseInt(x[2]),x[3],x[5])); }
        return out;
    }
    private static int[] readRanks(Path p,String id) throws Exception {
        final ArrayList<Integer> x=new ArrayList<Integer>(); final String[] ls=new String(Files.readAllBytes(p),UTF8).split("\\r?\\n",-1);
        for(String line:ls) { String s=line.trim(); if(s.length()==0 || s.startsWith("#")) continue; try{x.add(Integer.valueOf(s));}catch(NumberFormatException e){throw new IOException(id+": malformed rank "+s);}}
        int[] r=new int[x.size()]; for(int i=0;i<r.length;i++) r[i]=x.get(i); return r;
    }
    /** CellNetParser's constructor still mutates the shared ml.CellNet.DENSE global while
     * building a master net (it no longer matters for inference -- see InferenceNet -- but the
     * parser itself still writes it). Serialized on CellNet.class so two CONCURRENT bundle
     * loads/packs in this JVM cannot race each other's master construction. This does not, and
     * cannot, protect against some unrelated, uncoordinated caller elsewhere in the same JVM
     * also calling CellNetParser directly without this lock. */
    private static CellNet parseNet(byte[] bytes,String id) throws Exception {
        synchronized(CellNet.class) {
            final ArrayList<byte[]> lines=new ArrayList<byte[]>(); String[] ls=new String(bytes,UTF8).split("\\n",-1);
            for(String s:ls) { if(s.endsWith("\r")) s=s.substring(0,s.length()-1); lines.add(s.getBytes(UTF8)); }
            try { CellNet n=CellNetParser.loadFromLines(lines); if(n==null) throw new IOException(id+": null CellNet"); return n; }
            catch(Throwable e) { if(e instanceof IOException) throw (IOException)e; throw new IOException(id+": CellNet parse failure",e); }
        }
    }
    /** The exact header-region boundary CellNetParser.parseHeader() uses (ml/CellNetParser.java:93-163):
     * blank lines are skipped WITHIN the header, and scanning stops at the first line that is both
     * non-empty and does not start with "#" (the first cell/edge/weight data line) -- never
     * scanning past it. Every header helper below reads from this bounded slice, not the whole
     * net-bytes content, so a "#dense"/"#density"/etc.-shaped prefix could never accidentally match
     * text living in the data section. In practice this made no observable difference for any
     * currently well-formed net (data lines start with 'C'/'W'/'I'/'H' or a double-hash "##"
     * comment, never a bare single-hash header prefix) -- but it makes the boundary correct by
     * construction instead of correct by accident. Root's review, 2026-09-09: the prior versions of
     * these helpers scanned the entire file. */
    private static String[] headerLines(byte[] bytes) {
        final String[] all=new String(bytes,UTF8).split("\\r?\\n");
        int end=0;
        while(end<all.length) {
            final String l=all[end];
            if(l.length()==0) { end++; continue; }
            if(!l.startsWith("#")) break;
            end++;
        }
        return Arrays.copyOfRange(all,0,end);
    }
    private static int[] parseDims(byte[] bytes,String id) throws Exception {
        for(String line:headerLines(bytes)) if(line.startsWith("#dims")) {
            String s=line.substring(5).trim(); return parseInts(s,"\\s+",id+" dims");
        }
        throw new IOException(id+": missing #dims");
    }
    private static int parseSeed(byte[] bytes,String id) throws Exception {
        for(String l:headerLines(bytes)) if(l.startsWith("#seed")) return Integer.parseInt(l.substring(5).trim());
        throw new IOException(id+": missing #seed");
    }
    /** `ml.CellNet.DENSE` is a single GLOBAL static flag, set by CellNetParser's constructor
     * from whichever file it most recently parsed -- NOT a per-CellNet-instance property.
     * CellNet.feedForward() reads that global flag at CALL time to choose feedForwardDense()
     * vs feedForwardSparse(), so a bundle holding a MIX of dense and sparse subnets (this
     * candidate export is #sparse; every synthetic legacy subnet in this increment's own
     * fixtures is #dense) will silently score a later-scored net under the WRONG algorithm
     * once any other-format net is parsed afterward -- reproduced directly: isolated parse+score
     * of a real #sparse four-output candidate matched Torch to 5.96e-8; the same net scored
     * from inside a 176-subnet mixed bundle (where a #dense subnet is parsed last, during
     * pack()/validate()'s sequential loop) crashed with `AssertionError: 22, 1` inside
     * feedForwardDense() -- a shape mismatch from running dense summation on a sparsely-
     * structured Cell. Root cause confirmed by Yoimiya's source-read hypothesis, 2026-09-09;
     * empirically reproduced and isolated by Ady the same day (see
     * results/subnet_four_output_bundle_candidate_integration_20260909.md). A first fix (reset
     * the global flag immediately before every score()/scoreAll() call) was REJECTED: it raced
     * when different-format subnets were scored concurrently, and did not protect the direct
     * netForCurrentThread().feedForward() call path MagQCVectorMaker already uses. The accepted
     * replacement (see InferenceNet, above) removes the global write entirely -- each Subnet's
     * own dense/sparse-ness, captured here, is instead baked into a per-instance feedForward()
     * override, so there is no shared mutable state left to race or bypass.
     * <p>Defaults to TRUE when neither #dense nor #sparse appears in the header -- matching
     * CellNetParser.parseHeader()'s own field default (`boolean dense=true;`, ml/CellNetParser.java:334)
     * exactly, rather than rejecting a net the real production parser would load fine. The prior
     * version of this method threw on absence, which was STRICTER than the actual
     * producer/consumer contract -- a real, if so-far-latent (every net used this session always
     * carried an explicit marker), correctness gap: root's review, 2026-09-09, verified directly
     * against CellNetParser source, not assumed. */
    private static boolean parseDenseFlag(byte[] bytes,String id) throws Exception {
        boolean dense=true; // CellNetParser's own default when neither marker is present
        for(String l:headerLines(bytes)) {
            // Last matching line wins WITHIN the (now correctly bounded) header region -- matches
            // CellNetParser.parseHeader()'s own loop, which reassigns `dense` on every #dense/
            // #sparse line with no early exit for repeats (it DOES stop at the first non-header
            // line -- see headerLines()).
            if(l.startsWith("#dense")) dense=true;
            else if(l.startsWith("#sparse")) dense=false;
        }
        return dense;
    }
    /** density/density1/blockSize, matching CellNetParser's own defaults (1f/0f/1 when absent)
     * and header precedence exactly: `#density1` is checked BEFORE `#density` because
     * "#density1".startsWith("#density") is true -- checking `#density` first would silently
     * misread a `#density1` line as `#density`, which is exactly the kind of misinterpretation
     * this helper exists to rule out. Last matching line wins, same reasoning as
     * parseDenseFlag(). */
    private static float[] parseDensities(byte[] bytes,String id) throws Exception {
        float density=1f, density1=0f;
        for(String l:headerLines(bytes)) {
            if(l.startsWith("#density1")) density1=Float.parseFloat(l.substring(9).trim());
            else if(l.startsWith("#density")) density=Float.parseFloat(l.substring(8).trim());
        }
        return new float[]{density,density1};
    }
    private static int parseBlockSize(byte[] bytes,String id) throws Exception {
        int blockSize=1;
        for(String l:headerLines(bytes)) {
            if(l.startsWith("#blocksize")) blockSize=Integer.parseInt(l.substring(10).trim());
        }
        return blockSize;
    }

    private static Path resolve(Path root,String supplied) throws IOException {
        if("-".equals(supplied)) throw new IOException("cannot resolve '-'");
        Path p=Paths.get(supplied); if(Files.isRegularFile(p)) return p;
        Path q=root.resolve(p.getFileName().toString()); if(Files.isRegularFile(q)) return q;
        q=root.resolve("nets").resolve(p.getFileName().toString()); if(Files.isRegularFile(q)) return q;
        q=root.resolve("subsets").resolve(p.getFileName().toString()); if(Files.isRegularFile(q)) return q;
        throw new IOException("cannot resolve staged path "+supplied+" below "+root);
    }
    /** Accepts t/true or f/false (case-insensitive); absent means false; any other supplied
     * value is rejected outright rather than silently defaulting to false (a typo like "tru"
     * must not silently select the legacy loader). */
    private static boolean parseMultioutputFlag(String v) {
        if(v==null) return false;
        if("t".equalsIgnoreCase(v) || "true".equalsIgnoreCase(v)) return true;
        if("f".equalsIgnoreCase(v) || "false".equalsIgnoreCase(v)) return false;
        throw new IllegalArgumentException("multioutput="+v+": expected t/true or f/false");
    }
    private static Map<String,String> parseArgs(String[] a,int start) { Map<String,String> m=new LinkedHashMap<String,String>(); for(int i=start;i<a.length;i++){int e=a[i].indexOf('='); if(e<=0) throw new IllegalArgumentException("Expected key=value: "+a[i]); m.put(a[i].substring(0,e).toLowerCase(),a[i].substring(e+1));} return m; }
    private static Path requiredPath(Map<String,String> a,String k) { String v=a.get(k); if(v==null||v.length()==0) throw new IllegalArgumentException("Required: "+k+"="); return Paths.get(v); }
    private static Path pathOr(Map<String,String> a,String k,Path d) { return a.containsKey(k)?Paths.get(a.get(k)):d; }
    private static String req(Map<String,String> m,String k) throws IOException { String v=m.get(k); if(v==null) throw new IOException("missing field "+k); return v; }
    private static void require(Map<String,String> m,String k) throws IOException { req(m,k); }
    private static int parseInt(Map<String,String> m,String k) throws IOException { try{return Integer.parseInt(req(m,k));}catch(NumberFormatException e){throw new IOException("bad integer "+k,e);} }
    private static int[] parseInts(String s,String delim,String what) throws IOException { if(s.length()==0) throw new IOException("empty "+what); String[] x=s.split(delim,-1); int[] r=new int[x.length]; try{for(int i=0;i<r.length;i++)r[i]=Integer.parseInt(x[i].trim());}catch(NumberFormatException e){throw new IOException("bad "+what,e);} return r; }
    private static boolean safeKey(String s) { return s!=null && s.matches("[A-Za-z0-9_.-]+") && s.length()<=200; }
    private static byte[] decode(String s,String what) throws IOException { try{return Base64.getDecoder().decode(s);}catch(IllegalArgumentException e){throw new IOException("bad base64 "+what,e);} }
    /** Single-line, unwrapped base64(UTF-8) of the given strings joined by LF with no trailing
     * LF -- the exact encoding root specified for output_names/output_units, deliberately NOT
     * using writeWrapped()'s continuation-line mechanism (avoids extending continuation parsing
     * for a payload this short). */
    private static String b64Lines(String[] parts) {
        final StringBuilder sb=new StringBuilder();
        for(int i=0;i<parts.length;i++) { if(i>0) sb.append('\n'); sb.append(parts[i]); }
        return Base64.getEncoder().encodeToString(sb.toString().getBytes(UTF8));
    }
    private static String[] decodeLines(String b64,String what) throws IOException {
        return new String(decode(b64,what),UTF8).split("\n",-1);
    }
    private static String b64(byte[] x) { return Base64.getEncoder().encodeToString(x); }
    private static String sha256(byte[] x) { try{MessageDigest d=MessageDigest.getInstance(HASH_ALGORITHM); return hex(d.digest(x));}catch(Exception e){throw new RuntimeException(e);} }
    private static String hex(byte[] x){StringBuilder b=new StringBuilder(x.length*2);for(byte v:x)b.append(String.format("%02x",v&255));return b.toString();}
    private static String join(int[] x,String sep){StringBuilder b=new StringBuilder();for(int i=0;i<x.length;i++){if(i>0)b.append(sep);b.append(x[i]);}return b.toString();}
    private static String ranksToken(int[] x){return x==null?"-":(x.length==0?"empty":join(x,","));}
    private static String join(String[] x,String sep){StringBuilder b=new StringBuilder();for(int i=0;i<x.length;i++){if(i>0)b.append(sep);b.append(x[i]);}return b.toString();}

    private static final class ReleaseRow { final String id,type,subsetPath; final int expectedOutputs; ReleaseRow(String i,String t,String s,int o){id=i;type=t;subsetPath=s;expectedOutputs=o;} }
    private static final class AggRow { final String id,type,subsetPath,netPath; final int numObs,expectedInputs; AggRow(String i,String t,int o,int in,String s,String n){id=i;type=t;numObs=o;expectedInputs=in;subsetPath=s;netPath=n;} }
    private static final class SubnetText { final Map<String,String> fields; final String netB64; SubnetText(Map<String,String> f,String n){fields=f;netB64=n;} }
}
