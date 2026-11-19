package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaCenteredGatePool;
import prok.NcrnaCenteredGateReplay;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Six raw loci, including a reverse-strand third-member locus whose reference
 * name resembles the main member. Two by three tied choices must yield all six
 * clan combinations, with hand-derived integer sums and unchanged third calls.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGatePoolTest{
	/*--------------------------------------------------------------*/
	/*----------------          Entry Point         ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGatePoolTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown finite pooling argument");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null, "Finite pooling requires only out=fresh_directory");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fresh finite pooling root required");
			run(root);
			pp.outstream.println("NCRNA_CENTERED_GATE_POOL_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Finite Acceptance      ----------------*/
	/*--------------------------------------------------------------*/

	/** Includes a distant unchanged third-member start, so full and raw-near
	 * arithmetic differ while all six loci remain in whole-cohort outcomes. */
	private static void run(final File root){
		final ByteBuilder manifest=new ByteBuilder(NcrnaCenteredGatePool.MANIFEST_HEADER).nl();
		final ByteBuilder witnesses=new ByteBuilder(NcrnaCenteredGatePool.WITNESS_MANIFEST_HEADER).nl();
		for(final String panel:new String[]{"new1000", "old300"}){
			final String witness=path(root, panel+".witnesses.tsv");
			write(witness, rawWitnesses());
			witnesses.append(panel).tab().append(witness).append("\t2\t2\t2\n");
			for(int seed=58; seed<=60; seed++){
				for(int member=0; member<2; member++){
					final String prefix=path(root, panel+"_"+MEMBERS[member]+"_"+seed);
					final ByteBuilder rows=new ByteBuilder(NcrnaCenteredGateReplay.RESULT_HEADER).nl();
					for(final String policy:POLICIES){
						for(final boolean joint:new boolean[]{false, true}){
							for(int end=0; end<2; end++){
								final int variants=policy.equals("BEST") && (joint || end==0) ? member+2 : 1;
								for(int variant=0; variant<variants; variant++){
									rows.append(summary(panel, seed, member, policy, joint, end, variant)).nl();
								}
							}
						}
					}
					write(prefix+".summary.tsv", rows);
					final ByteBuilder unavailable=new ByteBuilder("member\tseed\tpolicy\tend\tside\tmonitor_target_moves\tqualification\n");
					if(member==1){unavailable.append("RF01482\t").append(seed).append("\tRATE_MARGIN0\t3prime\tHIGH_UNAVAILABLE\t3\tMONITOR_UNATTAINABLE_NOT_RETUNED_ON_OUTER\n");}
					write(prefix+".unavailable.tsv", unavailable);
					manifest.append(panel).tab().append(seed).tab().append(MEMBERS[member]).tab()
						.append(prefix).append(".summary.tsv\t").append(prefix).append(".unavailable.tsv\n");
				}
			}
		}
		write(path(root, "summaries.tsv"), manifest);
		write(path(root, "witnesses.tsv"), witnesses);
		new NcrnaCenteredGatePool(arguments(root, "summaries.tsv", "witnesses.tsv", "result")).process();
		int rows=0, jointBest=0, neverRows=0;
		final ObjectSet<String> identities=new ObjectSet<String>(String.class);
		final ByteBuilder identity=new ByteBuilder();
		try(Input in=new Input(path(root, "result/clan.tsv"))){
			final byte[] header=in.reader.nextLine();
			check(new String(header, java.nio.charset.StandardCharsets.US_ASCII).contains("refined_members_fixed"), "Unmeasured third-member physical status must not be represented as an assessed clan total");
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==53 && p.parseLong(13)==6, "The full clan includes all six loci");
				for(int i=0; i<p.terms(); i++){
					check(!p.termEquals("NA", i) && !p.termEquals("NaN", i) && !p.termEquals("Infinity", i),
						"Pooling must preserve conditional absence without emitting numeric placeholders");
				}
				for(int i=0; i<4; i++){
					final int side=i<2 ? 5+i : 7+i;
					check(p.length(side)==0 ? p.length(49+i)==0
						: p.termEquals(p.termEquals("RATE_MARGIN01", 2) ? "NEVER" : "THRESHOLD", 49+i),
						"Every pooled member/end retains its own active threshold mode");
				}
				identity.clear();
				for(int i=0; i<13; i++){
					identity.append(p.parseString(i)).tab();
				}
				check(identities.add(identity.toString()), "Every reported tied cutoff tuple has its own complete identity");
				final boolean five=p.termEquals("5prime", 4);
				if(p.termEquals("RATE_MARGIN01", 2)){
					check(p.parseLong(23)==0 && p.parseLong(24)==0 && p.parseLong(25)==6 && p.parseLong(26)==0
						&& p.parseLong(44)==p.parseLong(45) && p.parseLong(46)==p.parseLong(47),
						"Never-move tuples retain all six raw loci, errors and zero movements");
					neverRows++;
				}
				check(p.parseLong(28)==(five ? 5 : 6), "The distant third-member start remains in the full cohort but outside the raw-near subset");
				check(p.parseLong(44)==(five ? 183 : 25) && p.parseLong(46)==(five ? 14 : 25), "Distinct exact full and near square sums include the proper unchanged third-member contributions");
				if(p.termEquals("BEST", 2)){
					check(p.parseLong(45)==(five ? 178 : 20) && p.parseLong(47)==(five ? 9 : 20)
						&& p.parseLong(14)==4 && p.parseLong(23)==2 && p.parseLong(24)==0 && p.parseLong(25)==4 && p.parseLong(26)==2,
						"Only two refined-member ends improve; the two third-member ends remain tied and unchanged");
					if(p.termEquals("PAIRED_JOINT", 3)){
						check(p.parseLong(17)==3, "Both-exact clan count includes the raw third-member exact pair");
						check((p.parseFloat(7)==.1f || p.parseFloat(7)==.2f)
							&& (p.parseFloat(11)==.1f || p.parseFloat(11)==.2f || p.parseFloat(11)==.3f),
							"Both members' tied thresholds stay explicit rather than being collapsed into an anonymous pooled point");
						jointBest++;
					}else{check(p.length(17)==0, "Single-end tables do not fabricate joint measurements");}
				}
				rows++;
			}
		}
		check(rows==306 && jointBest==72 && neverRows==24, "Every panel/seed retains the complete two-by-three tied-cutoff product at both joint ends");
		int absent=0;
		try(Input in=new Input(path(root, "result/unavailable.tsv"))){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==8 && p.termEquals("RF01482", 1) && p.termEquals("HIGH_UNAVAILABLE", 5) && p.parseInt(6)==3,
					"Unavailable monitor targets remain visible with panel/member/seed identity");
				absent++;
			}
		}
		check(absent==6, "All six saved unavailable brackets survive reporting without invented matching cutoffs");
		legacy(root, manifest);
		final ByteBuilder duplicate=rawWitnesses();
		try(Input in=new Input(path(root, "new1000.witnesses.tsv"))){
			in.reader.nextLine();
			duplicate.append(in.reader.nextLine()).nl();
		}
		write(path(root, "duplicate.witnesses.tsv"), duplicate);
		write(path(root, "bad_witnesses.tsv"), new ByteBuilder(witnesses.toString().replace(path(root, "new1000.witnesses.tsv"), path(root, "duplicate.witnesses.tsv"))));
		expectFailure(arguments(root, "summaries.tsv", "bad_witnesses.tsv", "bad_duplicate"), "reference cannot be counted twice");
		final String original=path(root, "new1000_RF00174_58.summary.tsv");
		final ByteBuilder wrong=new ByteBuilder();
		try(Input in=new Input(original)){
			wrong.append(in.reader.nextLine()).nl();
			final LineParser1 p=new LineParser1('\t');
			boolean first=true;
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				if(first){
					p.set(line);
					for(int i=0; i<48; i++){
						if(i>0){wrong.tab();}
						if(i==41){wrong.append(p.parseLong(i)+1);}
						else if(i==17){wrong.append(Double.toString(Math.sqrt((p.parseLong(41)+1)/2.0)));}
						else if(i==19){wrong.append(Double.toString(Math.sqrt(p.parseLong(42)/2.0)-Math.sqrt((p.parseLong(41)+1)/2.0)));}
						else{wrong.append(p.parseString(i));}
					}
					wrong.nl();
					first=false;
				}else{wrong.append(line).nl();}
			}
		}
		write(path(root, "wrong.summary.tsv"), wrong);
		write(path(root, "bad_summaries.tsv"), new ByteBuilder(manifest.toString().replace(original, path(root, "wrong.summary.tsv"))));
		expectFailure(arguments(root, "bad_summaries.tsv", "witnesses.tsv", "bad_cohort"), "exact original raw cohort");
	}

	/** The complete old header permits old optional cells, normalized before any output. */
	private static void legacy(final File root, final ByteBuilder manifest){
		String legacyManifest=manifest.toString();
		final ByteBuilder text=new ByteBuilder();
		final LineParser1 p=new LineParser1('\t');
		for(final String panel:new String[]{"new1000", "old300"}){
			for(int seed=58; seed<=60; seed++){
				for(int member=0; member<2; member++){
					final String prefix=path(root, panel+"_"+MEMBERS[member]+"_"+seed);
					final String current=prefix+".summary.tsv", old=prefix+".legacy.tsv";
					legacyManifest=legacyManifest.replace(current, old);
					try(Input in=new Input(current); Output out=new Output(old)){
						in.reader.nextLine();
						out.writer.println(NcrnaCenteredGateReplay.LEGACY_RESULT_HEADER);
						for(byte[] row; (row=in.reader.nextLine())!=null;){
							p.set(row);
							check(p.terms()==48, "Legacy fixture conversion requires the complete new summary");
							text.clear();
							for(int i=0; i<46; i++){
								if(i>0){text.tab();}
								if((i==6 || i==7) && p.termEquals("NEVER", 40+i)){text.append("Infinity");}
								else{text.append(p.length(i)==0 ? "NA" : p.parseString(i));}
							}
							out.writer.print(text.nl());
						}
					}
				}
			}
		}
		write(path(root, "legacy_summaries.tsv"), new ByteBuilder(legacyManifest));
		new NcrnaCenteredGatePool(arguments(root, "legacy_summaries.tsv", "witnesses.tsv", "legacy_result")).process();
		try(Input current=new Input(path(root, "result/clan.tsv")); Input old=new Input(path(root, "legacy_result/clan.tsv"))){
			for(byte[] row; (row=current.reader.nextLine())!=null;){
				check(Arrays.equals(row, old.reader.nextLine()), "Legacy normalization cannot change counts, sums, modes or tied tuples");
			}
			check(old.reader.nextLine()==null, "Legacy pooling cannot add report rows");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Literal Evidence       ----------------*/
	/*--------------------------------------------------------------*/

	private static ByteBuilder rawWitnesses(){
		final ByteBuilder text=new ByteBuilder(NcrnaCenteredGatePool.WITNESS_HEADER).nl();
		for(int member=0; member<3; member++){
			for(int locus=0; locus<2; locus++){
				final int rs=1000+200*(member*2+locus), re=rs+119;
				final int d5=ERRORS[member][locus][0], d3=ERRORS[member][locus][1];
				final boolean minus=locus==0;
				text.append("input\tRF00174_ref_").append(member).append('_').append(locus).tab().append(MEMBERS[member])
					.append("\tcobalamin_").append(MEMBERS[member]).append("_model\tcontig\t")
					.append(minus ? rs-d3 : rs+d5).tab().append(minus ? re-d5 : re+d3).tab().append(rs).tab().append(re)
					.tab().append(minus ? '-' : '+').tab().append(d5).tab().append(d3).nl();
			}
		}
		return text;
	}

	private static String summary(final String panel, final int seed, final int member, final String policy,
			final boolean joint, final int end, final int variant){
		final String[] f=new String[48];
		Arrays.fill(f, "0");
		f[0]=MEMBERS[member];
		f[1]=Integer.toString(seed);
		f[2]=MEMBERS[member]+"_"+panel;
		f[3]=policy;
		Arrays.fill(f, 4, 8, "");
		f[46]=f[47]="";
		final boolean gate=policy.equals("BEST") || policy.startsWith("RATE_");
		if(gate){
			for(int e=0; e<2; e++){
				if(!joint && e!=end){continue;}
				f[4+e]=policy.equals("BEST") ? "ALL_TIED_BEST" : "EXACT";
				final boolean never=policy.equals("RATE_MARGIN01");
				f[6+e]=never ? "" : Float.toString(e==0 ? .1f*(variant+1) : .6f);
				f[46+e]=never ? "NEVER" : "THRESHOLD";
			}
		}
		f[8]=joint ? "PAIRED_JOINT" : "SINGLE_OTHER_RAW";
		f[9]=end==0 ? "5prime" : "3prime";
		final long[] n=new long[46];
		long before=0, after=0;
		for(int locus=0; locus<2; locus++){
			final int raw=ERRORS[member][locus][end];
			final int otherRaw=ERRORS[member][locus][1-end];
			final int value=policy.equals("BEST") && member==end ? 0 : raw;
			final int other=joint && policy.equals("BEST") && member==1-end ? 0 : otherRaw;
			final int a=Math.abs(value), b=Math.abs(raw), z=Math.abs(other);
			n[10]++;
			if(a==0){n[11]++;}
			if(a<=3){n[12]++;}
			if(a<=10){n[13]++;}
			if(a==0 && z==0){n[14]++;}
			if(a<=3 && z<=3){n[15]++;}
			if(a<=10 && z<=10){n[16]++;}
			if(a<b){
				n[20]++;
			}else if(a>b){
				n[21]++;
			}else{n[22]++;}
			if(value!=raw){n[23]++;}
			if(b==0 && a>0){n[24]++;}
			n[25]++;
			if(a>b){n[28]++;}
			before+=(long)raw*raw;
			after+=(long)value*value;
		}
		for(int i=10; i<=39; i++){
			f[i]=Long.toString(n[i]);
		}
		if(!joint){Arrays.fill(f, 14, 17, "");}
		f[17]=Double.toString(Math.sqrt(before/2.0));
		f[18]=Double.toString(Math.sqrt(after/2.0));
		f[19]=Double.toString(Math.sqrt(after/2.0)-Math.sqrt(before/2.0));
		f[26]=f[17];
		f[27]=f[18];
		f[40]=n[20]+n[21]==0 ? "" : Double.toString((double)n[20]/(n[20]+n[21]));
		f[41]=Long.toString(before);
		f[42]=Long.toString(after);
		f[43]=f[41];
		f[44]=f[42];
		f[45]="CONDITIONAL_REAL_CALLS_NOT_CALLER_OR_NESTED_MONITOR";
		return String.join("\t", f);
	}
	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private static void expectFailure(final String[] args, final String message){
		boolean failed=false;
		try{
			new NcrnaCenteredGatePool(args).process();
		}catch(IllegalArgumentException expected){
			failed=expected.getMessage().contains(message);
		}
		check(failed, "Malformed saved evidence must fail specifically: "+message);
	}
	private static String[] arguments(final File root, final String a, final String b, final String out){
		return new String[]{"in="+path(root, a), "in2="+path(root, b), "out="+path(root, out), "t=1"};
	}
	private static String path(final File root, final String name){return new File(root, name).getAbsolutePath();}
	private static void write(final String path, final ByteBuilder content){
		try(Output out=new Output(path)){
			out.writer.print(content);
		}
	}
	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
	/*--------------------------------------------------------------*/
	/*----------------        Fixture Constants     ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final int[][][] ERRORS={{{2,0},{-1,4}},{{3,-2},{0,1}},{{-13,2},{0,0}}};
	private static final String[] POLICIES={"RAW", "CENTER_ONLY", "CENTER_MARGIN0", "CENTER_MARGIN01", "ALL_WHOLE_PROPOSALS", "IGNORE_HEAD_DIAGNOSTIC", "BEST", "RATE_MARGIN0", "RATE_MARGIN01"};
}
