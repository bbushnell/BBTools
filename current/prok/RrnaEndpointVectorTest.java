package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Random;
import idaligner.AlignmentStats;

/** Synthetic nonzero-anchor/end-coordinate and measured-jitter fixtures.
 * @author Ganyu
 */
public final class RrnaEndpointVectorTest {

	public static void main(String[] args) throws Exception{
		final Path root=args.length==0?Files.createTempDirectory("rrna_endpoint_vector_"):Files.createDirectory(java.nio.file.Paths.get(args[0]));
		final byte[] bases=new byte[180];Arrays.fill(bases,(byte)'A');
		System.arraycopy("CGTACG".getBytes(StandardCharsets.US_ASCII),0,bases,47,6);System.arraycopy("TGCATG".getBytes(StandardCharsets.US_ASCII),0,bases,125,6);
		final RrnaPositionalKmerTable.Table five=new RrnaPositionalKmerTable.Table("test","5prime",6,-3,.5),three=new RrnaPositionalKmerTable.Table("test","3prime",6,-4,.5);
		for(int i=0;i<100;i++){five.add(bases,50);three.add(bases,129);}five.finish();three.finish();
		final float[] scratch=new float[25],a=new float[28],b=new float[28];
		check(RrnaEndpointVector.fill(five,bases,53,125,80,.4f,.75f,scratch,a),"Complete planted start window must be usable");
		check(RrnaPositionalKmerTableTest.best(Arrays.copyOf(a,25))==9,"Start signal at47 uses raw53 plus anchor-3; peak class9 must represent true end50, not anchor47");
		check(RrnaEndpointVector.endpointFromClass(53,9)==50,"Start label must recover planted inclusive gene end with nonzero anchor");
		check(RrnaEndpointVector.fill(three,bases,53,125,80,.4f,.75f,scratch,b),"Complete planted stop window must be usable");
		check(RrnaPositionalKmerTableTest.best(Arrays.copyOf(b,25))==16,"Stop signal at125 uses raw125 plus anchor-4; peak class16 must represent true end129");
		check(RrnaEndpointVector.endpointFromClass(125,16)==129,"Stop label must recover last gene base, not first base of anchor k-mer");
		// The source-contig example GGCCAAAAAA has4 GC bases/10. Candidate length is125-53+1=73; reference length80.
		check(a[25]==73/80f&&a[26]==4/10f&&a[27]==.75f,"Extras must be73/80 length ratio, full-contig4/10 GC and supplied candidate identity in that order");
		check(b[25]==a[25]&&b[26]==a[26]&&b[27]==a[27],"Two ends of one raw candidate share all three extras");
		final byte[] consensus=Arrays.copyOfRange(bases,50,130);final AlignmentStats stats=new AlignmentStats(true);
		check(RrnaEndpointVector.candidateIdentity(bases,50,129,consensus,stats)==1f,"Self-alignment must return identity1 from the exact inclusive candidate slice");
		final byte[] lower=new String(bases,StandardCharsets.US_ASCII).toLowerCase(java.util.Locale.ROOT).getBytes(StandardCharsets.US_ASCII);
		final byte[] mixed=consensus.clone();for(int i=0;i<mixed.length;i+=2){mixed[i]=(byte)(mixed[i]+32);}final byte[] lowerBefore=lower.clone(),mixedBefore=mixed.clone();
		check(RrnaEndpointVector.candidateIdentity(lower,50,129,mixed,stats)==1f,"Lowercase candidate and mixed-case reference must equal uppercase identity");
		check(Arrays.equals(lower,lowerBefore)&&Arrays.equals(mixed,mixedBefore),"Identity normalization must not mutate source contig or shared reference");
		stats.rStart=100000;stats.rStop=100100;stats.matches=999;stats.matchString="dirty".getBytes(StandardCharsets.US_ASCII);stats.doTrace=false;
		final AlignmentStats fresh=new AlignmentStats(true);final float expectedDirty=RrnaEndpointVector.candidateIdentity(bases,53,125,consensus,fresh);
		check(RrnaEndpointVector.candidateIdentity(lower,53,125,mixed,stats)==expectedDirty&&stats.rStart==fresh.rStart&&stats.matches==fresh.matches&&Arrays.equals(stats.matchString,fresh.matchString),"Dirty reused stats must match a fresh alignment including traceback, independent of candidate order");
		check(!RrnaEndpointVector.fill(five,bases,0,70,80,.4f,1f,scratch,b)&&Float.isNaN(b[0]),"Clipped scores cannot silently become neutral training features");
		final byte[] ambiguous=bases.clone();ambiguous[47]='N';check(!RrnaEndpointVector.fill(five,ambiguous,53,125,80,.4f,1f,scratch,b),"Ambiguous anchor context requires explicit exclusion, not fabricated scores");
		check(RrnaEndpointVector.integerField("r\tlflank=50 rflank=50 gccontig=0.4","lflank")==50,"Header parser must support whitespace separators");
		check(RrnaEndpointVector.floatField("r gccontig=0.4","gccontig")==.4f,"GC field must parse without a trailing space");
		testHistogram(root);
		final String first=root.resolve("first.tsv.gz").toString(),last=root.resolve("last.tsv.gz").toString();five.write(first);three.write(last);
		final Path input=root.resolve("loci.fa"),ref=root.resolve("refs.fa"),hist=root.resolve("hist.tsv");
		Files.write(input,(">a model=test lflank=50 rflank=50 gccontig=0.4\n"+new String(bases,StandardCharsets.US_ASCII)+"\n>b model=test lflank=50 rflank=50 gccontig=0.4\n"+new String(bases,StandardCharsets.US_ASCII)+'\n').getBytes(StandardCharsets.US_ASCII));
		Files.write(ref,(">test\n"+new String(consensus,StandardCharsets.US_ASCII)+'\n').getBytes(StandardCharsets.US_ASCII));
		Files.write(hist,"model\tend\tsigned_error\tcount\ntest\t5prime\t3\t20\ntest\t3prime\t-4\t20\n".getBytes(StandardCharsets.US_ASCII));
		final String out=root.resolve("vectors").toString();final String[] cli={"in="+input,"model=test","consensus="+ref,"starttable="+first,"stoptable="+last,"histogram="+hist,"outprefix="+out,"expected=2","samples=3","seed=74003"};
		RrnaEndpointVector.main(cli);
		final float id=RrnaEndpointVector.candidateIdentity(bases,53,125,consensus,stats);
		RrnaEndpointVector.fill(five,bases,53,125,80,.4f,id,scratch,a);RrnaEndpointVector.fill(three,bases,53,125,80,.4f,id,scratch,b);
		checkRows(out+".5prime.tsv.gz",a,9,53,50);checkRows(out+".3prime.tsv.gz",b,16,125,129);
		final String audit=RrnaEndpointTestSupport.read(out+".audit.tsv.gz");check(audit.split("\n").length==13&&audit.contains("\t50\t50\t9\tKEPT\t0"),"Audit must bind the same start truth50, anchorcenter50, label9 and first training row");
		final String summary=RrnaEndpointTestSupport.read(out+".stats.tsv");check(summary.contains("test\t5prime\t2\t6\t6\t0\t0\t20")&&summary.contains("test\t3prime\t2\t6\t6\t0\t0\t20"),"Expected2loci*3trials must be conserved per end");
		boolean failed=false;try{RrnaEndpointVector.main(cli);}catch(IllegalArgumentException e){failed=true;}check(failed,"A rerun must not overwrite the evidence from a completed dataset");
		System.out.println("PASS RrnaEndpointVectorTest nonzero_anchors both_ends candidate_extras case_normalization input_nonmutation dirty_stats histogram_clamping seeded_determinism missing_histogram missing_scores cli_api_parity labels counts fresh_outputs");
	}
	static void checkRows(String path,float[] features,int label,int raw,int truth) throws Exception{
		final String text=RrnaEndpointTestSupport.read(path);check(text.startsWith("#dims\t28\t25\n")&&text.contains("length_ratio\tcontig_gc\tidentity\tend_offset_-12"),"File header must declare trainer width, exact extras order and output offset order");int rows=0;
		for(String line:text.split("\n")){if(line.startsWith("#")){continue;}final String[] f=line.split("\t");check(f.length==53,"Every trainer row has28 inputs and25 targets");
			for(int j=0;j<28;j++){check(Math.abs(Float.parseFloat(f[j])-features[j])<1e-6,"CLI row must agree with shared feature API at input"+j);}
			int positives=0,best=-1;for(int j=0;j<25;j++){final int value=Integer.parseInt(f[28+j]);check(value==0||value==1,"One-hot labels are binary");if(value==1){positives++;best=j;}}
			check(positives==1&&best==label&&RrnaEndpointVector.endpointFromClass(raw,best)==truth,"Argmax of one-hot label must recover the independently planted true end");rows++;}
		check(rows==6,"Two loci with three trials each must yield six complete rows");
	}
	static void testHistogram(Path root) throws Exception{
		final RrnaEndpointVector.Distribution d=new RrnaEndpointVector.Distribution();d.add(-30,1);d.add(-12,2);d.add(0,3);d.add(12,4);d.add(50,5);
		check(d.total==15&&d.low==1&&d.high==5&&d.weights[0]==3&&d.weights[12]==3&&d.weights[24]==9,"Clipping transfers tail mass to±12 without dropping or rescaling counts");
		final Random x=new Random(71),y=new Random(71);for(int i=0;i<100;i++){final int a=d.draw(x);check(a==d.draw(y)&&(a==-12||a==0||a==12),"Seeded draws reproduce and never invent a bin absent from the measured clipped distribution");}
		final Path bad=root.resolve("missing_hist.tsv");Files.write(bad,"model\tend\tsigned_error\tcount\npooled\t5prime\t0\t1\npooled\t3prime\t0\t1\n".getBytes(StandardCharsets.US_ASCII));
		boolean failed=false;try{RrnaEndpointVector.distributions(bad.toString(),"test");}catch(IllegalArgumentException e){failed=true;}check(failed,"Pooled rows cannot stand in for missing measured model-specific histograms");
		failed=false;try{d.add(50,1);}catch(IllegalArgumentException e){failed=true;}check(failed,"Duplicate signed bins must not silently double sample weight");
	}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
