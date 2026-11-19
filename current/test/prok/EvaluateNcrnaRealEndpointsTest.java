package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import dna.AminoAcid;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prok.EvaluateNcrnaRealEndpoints;
import prok.FamilyResourceResolver;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFixedEndpointTrials;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO;
import shared.Shared;
import structures.ByteBuilder;

/** Small source-qualified real-schema fixtures using the declared existing payloads.
 * No fits, new resource tables, library copies or genomic population measurements.
 * Literal signed coordinates and malformed joins test the new adapter, while
 * saved-feature and score equality test strand and source identity preservation.
 * @author Brian, Ganyu
 */
public final class EvaluateNcrnaRealEndpointsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, EvaluateNcrnaRealEndpointsTest.class, false);
		try{
			final File root=new File(value(pp.args, "out"));
			check(root.mkdir(), "Finite fixture output directory must be fresh");
			pureControls();
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
			options.resources.root=value(pp.args, "familyroot");
			options.resources.consensus=value(pp.args, "ref");
			options.resources.models=value(pp.args, "hbm");
			options.resources.seeds=value(pp.args, "seeds");
			options.boundaryOverride=false;
			options.endpoints=false;
			final FamilyResourceResolver resolver=new FamilyResourceResolver();
			final String familyName=option(pp.args, "family", "cobalamin_rf00174");
			resolver.bindRoot(familyName, options.resources.root);
			final NcrnaFamily family=NcrnaFamilyLoader.load(familyName, options, resolver);
			check(family.library.length>1 && (!option(pp.args, "mode", "evaluate").equals("evaluate")
				|| family.library[0].length>=112 && family.library[1].length>=109),
				"Finite controls need two real paired models; inference controls also require ordinary-length cores");
			fixture(root, family, pp.args);
			pp.outstream.println("NCRNA_REAL_ENDPOINT_ASSAY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void pureControls(){
		final ByteBuilder maximum=new ByteBuilder("left\t");
		NcrnaTextIO.appendMaximumLength(maximum, Integer.MAX_VALUE);
		check(maximum.append("\tright").toString().equals("left\t\tright"), "An unbounded maximum has no written number");
		maximum.clear();
		NcrnaTextIO.appendMaximumLength(maximum, 200);
		check(maximum.toString().equals("200"), "A finite maximum keeps its exact decimal value");
		check(maximum("", false)==Integer.MAX_VALUE && maximum("200", false)==200
			&& maximum("2147483647", true)==Integer.MAX_VALUE, "Current absence and the declared legacy control preserve in-memory geometry");
		for(final String bad:new String[]{"2147483647", "2147483648", "4294967496", "-1", "0", "NA", "1x", "999999999999999999999999"}){
			fails(()->maximum(bad, false));
		}
		fails(()->maximum("", true));
		check(EvaluateNcrnaRealEndpoints.senseStart(11, 50, 100, false)==10
			&& EvaluateNcrnaRealEndpoints.senseStop(11, 50, 100, false)==49
			&& EvaluateNcrnaRealEndpoints.senseStart(11, 50, 100, true)==50
			&& EvaluateNcrnaRealEndpoints.senseStop(11, 50, 100, true)==89, "Both inclusive genomic frames have literal expected coordinates");
		check(EvaluateNcrnaRealEndpoints.fixedMask("RNA,model:m n=7;partial=10", false)==1
			&& EvaluateNcrnaRealEndpoints.fixedMask("RNA,model:m n=7;partial=10", true)==2
			&& EvaluateNcrnaRealEndpoints.fixedMask("RNA;partial=01", true)==1
			&& EvaluateNcrnaRealEndpoints.fixedMask("RNA;partial=11", true)==3
			&& EvaluateNcrnaRealEndpoints.fixedMask("RNA,model:partial=10", false)==0, "Only exact partial keys fix ends; minus swaps the two genomic bits");
		check("m n=7".equals(EvaluateNcrnaRealEndpoints.attribute("RNA,model:m n=7,family:cobalamin_rf00174", "model")), "Full model names retain abundance metadata");
		fails(()->EvaluateNcrnaRealEndpoints.fixedMask("RNA;partial=10;partial=01", false));
		fails(()->EvaluateNcrnaRealEndpoints.fixedMask("RNA;partial=00", false));
		check(EvaluateNcrnaRealEndpoints.geometry(10, 9, 100, 20, 60)==1
			&& EvaluateNcrnaRealEndpoints.geometry(-1, 30, 100, 20, 60)==2
			&& EvaluateNcrnaRealEndpoints.geometry(10, 28, 100, 20, 60)==3
			&& EvaluateNcrnaRealEndpoints.geometry(10, 70, 100, 20, 60)==4
			&& EvaluateNcrnaRealEndpoints.geometry(10, 29, 100, 20, 60)==0, "Joint holds distinguish cross, contig bounds and both inclusive length limits");
		final float[] score=new float[22];
		score[7]=1;
		check(Arrays.equals(EvaluateNcrnaRealEndpoints.proposal(30, 160, new float[][]{score, null}, 0), new int[]{27, 160}), "An unavailable end holds while the other uses its signed proposal");
		score[10]=.95f;
		check(Arrays.equals(EvaluateNcrnaRealEndpoints.proposal(30, 160, new float[][]{score, null}, .1f), new int[]{30, 160}), "The second margin reuses scores and rejects the insufficient gap");
		score[21]=2;
		check(EvaluateNcrnaRealEndpoints.proposal(30, 160, new float[][]{score, null}, 0)[0]==30, "The out-of-range head vetoes rather than becoming a coordinate");
		final float[][] inward={new float[22], new float[22]};
		inward[0][20]=1;
		inward[1][0]=1;
		final int[] pair=EvaluateNcrnaRealEndpoints.proposal(30, 149, inward, 0);
		check(pair[0]==40 && pair[1]==139
			&& EvaluateNcrnaRealEndpoints.geometry(pair[0], pair[1], 200, 109, 200)==3
			&& EvaluateNcrnaRealEndpoints.geometry(pair[0], 149, 200, 109, 200)==0
			&& EvaluateNcrnaRealEndpoints.geometry(30, pair[1], 200, 109, 200)==0,
			"Each ten-base inward move is valid alone, but their paired span is too short; scopes cannot be silently interchanged");
	}

	private static void fixture(final File root, final NcrnaFamily family, final String[] original){
		final String familyName=option(original, "family", "cobalamin_rf00174");
		final String member=familyName.equals("cobalamin_rf00174") ? "RF00174" : "RF01482";
		final boolean gate=option(original, "mode", "evaluate").equals("gatevectors");
		final boolean geometry=gate || option(original, "mode", "evaluate").equals("geometry");
		final String fastaA=path(root, "a.fa"), fastaB=path(root, "b.fa"), list=path(root, "inputs.txt");
		final String roster=path(root, "sources.tsv"), gff=path(root, "calls.gff"), witnesses=path(root, "witnesses.tsv");
		final String offsets=path(root, "offsets.tsv");
		final ByteBuilder[] fasta={new ByteBuilder(), new ByteBuilder()};
		final ByteBuilder[] gffs={new ByteBuilder("##input_ordinal0\t0\n"), new ByteBuilder("##input_ordinal0\t1\n")};
		final ByteBuilder rows=new ByteBuilder("input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n");
		for(int i=0; i<7; i++){
			final int source=i<5 ? 0 : 1;
			final int model=i==6 ? 1 : 0;
			final boolean minus=i==1;
			final String name=i==0 || i==5 ? "shared" : "contig"+i;
			final int core=family.library[model].length;
			final byte[] bases=new byte[core+64];
			for(int j=0; j<bases.length; j++){bases[j]=(byte)(i==5 ? 'G' : "ACGT".charAt(j&3));}
			System.arraycopy(family.library[model], 0, bases, 32, core);
			int start=i==2 ? 0 : i==4 ? 2 : 32;
			int stop=start+core-1;
			int truthStart=i==2 || i==4 || i==6 ? start : start-5;
			int truthStop=i==4 || i==6 ? stop : stop-3;
			if(i==3){bases[start]='N';}
			if(gate && i==3){
				bases[start-5]='N';
				bases[stop-3]='N';
			}
			if(minus){AminoAcid.reverseComplementBasesInPlace(bases);}
			final int low=minus ? bases.length-stop : start+1;
			final int high=minus ? bases.length-start : stop+1;
			final int truthLow=minus ? bases.length-truthStop : truthStart+1;
			final int truthHigh=minus ? bases.length-truthStart : truthStop+1;
			fasta[source].append('>').append(name).nl().append(bases).nl();
			rows.append(source==0 ? "A" : "B").tab().append("ref").append(i).tab().append(member).tab()
				.append(family.modelNames[model]).tab().append(name).tab().append(low).tab().append(high).tab().append(truthLow).tab().append(truthHigh)
				.tab().append(minus ? "-" : "+").tab().append(start-truthStart).tab().append(stop-truthStop).nl();
			gffs[source].append(name).append("\tBBTools\tRNA\t").append(low).tab().append(high).append("\t1\t")
				.append(minus ? "-" : "+").append("\t0\tRNA,model:").append(family.modelNames[model]).append(",family:").append(familyName);
			if(i==2){gffs[source].append(";partial=10");}
			gffs[source].nl();
		}
		write(fastaA, fasta[0]);
		write(fastaB, fasta[1]);
		write(list, new ByteBuilder().append(fastaA).nl().append(fastaB).nl());
		write(roster, new ByteBuilder("input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n0\tassemblyA\tA\t")
			.append(fastaA).append("\tunusedA\n1\tassemblyB\tB\t").append(fastaB).append("\tunusedB\n"));
		write(gff, new ByteBuilder("##gff-version 3\n").append(gffs[0]).append(gffs[1]));
		write(witnesses, rows);
		final ByteBuilder offsetRows=new ByteBuilder(NcrnaFixedEndpointTrials.OFFSET_HEADER).nl();
		for(int model=0; model<family.modelNames.length; model++){
			offsetRows.append(member).tab().append(family.modelNames[model].split(" ", 2)[0])
				.append(model==0 ? "\t10\t10\t5\t3\n" : "\t0\t10\t0\t0\n");
		}
		write(offsets, offsetRows);
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(original));
		if(gate){
			args.clear();
			for(final String key:new String[]{"familyroot", "ref", "hbm", "seeds"}){args.add(key+"="+value(original, key));}
			args.add("family="+familyName);
			args.add("mode=geometry");
		}
		args.addAll(Arrays.asList("in="+witnesses, "gff="+gff, "sources="+roster, "sourceformat=native", "fastalist="+list,
			"inputs=2", "witnesses=7", "expected=7", "panel=finite", "out="+path(root, "result")));
		if(!geometry){args.add("modeloffsets="+offsets);}
		new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process();
		if(geometry){
			geometryRows(root, family, member);
			if(gate){
				gateFixture(root, original, familyName, offsets, roster, list);
				return;
			}
			final String legacy=path(root, "geometry_legacy.tsv");
			final String audit=path(root, "geometry_assembly_audit.tsv");
			write(legacy, new ByteBuilder("index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n")
				.append("11\tB\tGCF_B.1\t2\tbacteria\tunknown\tunusedB\t").append(fastaB).append("\t0\t0\t0\n")
				.append("29\tA\tGCF_A.1\t1\tbacteria\tunknown\tunusedA\t").append(fastaA).append("\t0\t0\t0\n"));
			write(audit, new ByteBuilder(prok.NcrnaEndpointFoldInputs.AUDIT_LABEL).nl()
				.append(prok.NcrnaEndpointFoldInputs.AUDIT_HEADER).nl()
				.append("/panel\t1\t").append(fastaA).append("\tassemblyA\tEXCLUDED_ASSEMBLY\n")
				.append("/panel\t2\t").append(fastaB).append("\tassemblyB\tEXCLUDED_ASSEMBLY\n"));
			args.add("sources="+legacy);
			args.add("sourceformat=legacy");
			args.add("out="+path(root, "geometry_legacy_result"));
			fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])));
			args.add("assemblyaudit="+audit);
			new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process();
			equalRows(path(root, "result.geometry.tsv"), path(root, "geometry_legacy_result.geometry.tsv"));
			final String badGff=path(root, "bad.gff");
			write(badGff, new ByteBuilder("##gff-version 3\n").append(gffs[0]).append(gffs[0]).append(gffs[1]));
			args.add("gff="+badGff);
			args.add("out="+path(root, "bad_duplicate"));
			fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process());
			check(!new File(path(root, "bad_duplicate.geometry.tsv")).exists(), "Ambiguous source provenance fails before geometry output");
			return;
		}
		final LineParser1 p=new LineParser1('\t');
		final String[] seedLabels=value(original, "seedlabels").split(",");
		final String[] newNets=value(original, "newnets").split(",");
		int expectedSummaryRows=0;
		for(int seed=0; seed<seedLabels.length; seed++){
			final boolean five=!newNets[2*seed].equals("none");
			final boolean three=!newNets[2*seed+1].equals("none");
			expectedSummaryRows+=8*((five ? 1 : 0)+(three ? 1 : 0))+(five && three ? 16 : 0);
		}
		int summaryRows=0;
		try(Input input=new Input(path(root, "result.summary.tsv"))){
			input.reader.nextLine();
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				check(p.terms()==44 && p.parseInt(5)==7 && p.parseInt(19)+p.parseInt(20)+p.parseInt(21)==7,
					"All metric rows preserve the seven calls and paired denominator");
				for(int field=0; field<p.terms(); field++){
					check(!p.termEquals("NA", field) && !p.termEquals("NaN", field) && !p.termEquals("Infinity", field),
						"Nonexistent statistics are omitted and measured assay values must remain finite");
				}
				check(p.parseLong(19)+p.parseLong(20)==0 ? p.length(22)==0
					: Math.abs(p.parseDouble(22)-(double)p.parseLong(19)/(p.parseLong(19)+p.parseLong(20)))<1e-7,
					"Movement-fraction presence and value follow the unchanged paired counts");
				int seed=0;
				while(seed<seedLabels.length && !p.termEquals(seedLabels[seed], 1)){seed++;}
				check(seed<seedLabels.length, "Summary seeds retain their declared identities");
				final int end=p.termEquals("5prime", 4) ? 0 : 1;
				check(!newNets[2*seed+end].equals("none"), "Unavailable centered models have no fabricated comparison rows");
				if(p.termEquals("SINGLE_END_OTHER_AT_BASE", 42)){
					check(p.length(11)==0 && p.length(12)==0, "Single-end comparisons cannot claim both-end accuracy");
				}else{
					check(p.termEquals("PAIRED_JOINT", 42) && !newNets[2*seed].equals("none") && !newNets[2*seed+1].equals("none"),
						"Only two completed models make a paired comparison");
				}
				check(p.parseInt(37)==1 && p.parseInt(38)==(p.termEquals("5prime", 4) ? 1 : 0) && p.parseInt(41)==1,
					"Exactly one invalid center, one fixed start and one unsupported model remain explicit");
				if(p.termEquals("center_only", 3)){
					check(p.parseInt(6)==7, "Literal training offsets make every fixture center exact; held and unsupported cases were already exact");
				}
				summaryRows++;
			}
		}
		check(summaryRows==expectedSummaryRows, "Every available end and complete pair has its declared four arms and two margins");
		try(Input input=new Input(path(root, "result.trials.tsv"))){
			input.reader.nextLine();
			for(byte[] row=input.reader.nextLine(); row!=null && row[0]!='#'; row=input.reader.nextLine()){
				p.set(row);
				check(p.terms()==28, "Trial provenance names the comparison scope and target end");
				if(!p.termEquals("SINGLE_END_OTHER_AT_BASE", 26)){continue;}
				final boolean centered=p.termEquals("center_only", 18) || p.termEquals("centered_nn", 18);
				final boolean validCenter=p.termEquals("VALID", 15);
				final int other=p.termEquals("5prime", 27) ? 1 : 0;
				final long base=p.parseLong((centered && validCenter ? 13 : 9)+other);
				check(p.parseLong(19+other)==base && p.parseLong(22+other)==base,
					"Every one-end proposal and final pair keeps the other coordinate at that arm's unmoved base");
			}
		}
		final String[] plus=new String[4], reverse=new String[4], otherSource=new String[4];
		int firstGc=0, secondGc=0;
		boolean missing=false, fixed=false, invalid=false;
		try(Input input=new Input(path(root, "result.features.tsv"))){
			input.reader.nextLine();
			for(byte[] row=input.reader.nextLine(); row!=null && row[0]!='#'; row=input.reader.nextLine()){
				p.set(row);
				final int index=(p.termEquals("raw", 1) ? 0 : 2)+(p.termEquals("5prime", 2) ? 0 : 1);
				if(p.parseInt(0)==1){plus[index]=new String(row, java.nio.charset.StandardCharsets.UTF_8).substring(2);}
				if(p.parseInt(0)==2){reverse[index]=new String(row, java.nio.charset.StandardCharsets.UTF_8).substring(2);}
				if(p.parseInt(0)==6){otherSource[index]=new String(row, java.nio.charset.StandardCharsets.UTF_8).substring(2);}
				if(index==0 && p.parseInt(0)==1){firstGc=p.parseInt(68);}
				if(index==0 && p.parseInt(0)==6){secondGc=p.parseInt(68);}
				missing|=p.parseInt(0)==4 && p.termEquals("MISSING_FEATURE_HOLD", 3);
				fixed|=p.parseInt(0)==3 && p.termEquals("FIXED", 3);
				invalid|=p.parseInt(0)==5 && p.termEquals("INVALID_CENTER_RAW_HOLD", 3);
			}
		}
		check(Arrays.equals(plus, reverse) && missing && fixed && invalid, "Plus/minus equivalents have identical feature bits, with missing/fixed/invalid-center states distinguished");
		check(!Arrays.equals(plus, otherSource) && firstGc!=secondGc,
			"Equal contig names and coordinates in two sources must retain their deliberately different genomic sequence/GC");
		final HashMap<String, String> plusScores=new HashMap<String, String>();
		int reverseScores=0;
		try(Input input=new Input(path(root, "result.scores.tsv"))){
			input.reader.nextLine();
			for(byte[] row=input.reader.nextLine(); row!=null && row[0]!='#'; row=input.reader.nextLine()){
				p.set(row);
				if(p.parseInt(0)!=1 && p.parseInt(0)!=2){continue;}
				check(p.terms()==27 && p.termEquals("AVAILABLE", 4) || p.terms()==5 && p.termEquals("MODEL_UNAVAILABLE", 4),
					"Equivalent interior calls produce every available score or an explicit missing-model record");
				final String key=p.parseString(1)+'\t'+p.parseString(2)+'\t'+p.parseString(3);
				final String content=new String(row, java.nio.charset.StandardCharsets.UTF_8).substring(2);
				if(p.parseInt(0)==1){check(plusScores.put(key, content)==null, "Each seed/arm/end is inferred once");}
				else{
					check(content.equals(plusScores.get(key)), "Equivalent plus/minus calls have identical score bits under the same seed/arm/end identity");
					reverseScores++;
				}
			}
		}
		check(reverseScores==4*value(original, "seedlabels").split(",").length, "Every seed/arm/end score pair was compared");
		final String legacy=path(root, "legacy.tsv");
		write(legacy, new ByteBuilder("index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n")
			.append("11\tB\tassemblyB\t2\tbacteria\tunknown\tunusedB\t").append(fastaB).append("\t0\t0\t0\n")
			.append("29\tA\tassemblyA\t1\tbacteria\tunknown\tunusedA\t").append(fastaA).append("\t0\t0\t0\n"));
		args.add("sources="+legacy);
		args.add("sourceformat=legacy");
		args.add("out="+path(root, "legacy_result"));
		new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process();
		for(final String suffix:new String[]{"features.tsv", "scores.tsv", "trials.tsv", "summary.tsv"}){
			equalRows(path(root, "result."+suffix), path(root, "legacy_result."+suffix));
		}
		final String badGff=path(root, "bad.gff");
		write(badGff, new ByteBuilder("##gff-version 3\n").append(gffs[0]).append(gffs[0]).append(gffs[1]));
		args.add("gff="+badGff);
		args.add("out="+path(root, "bad_duplicate"));
		fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process());
		check(!new File(path(root, "bad_duplicate.features.tsv")).exists(), "Ambiguous source provenance fails before output");
	}

	/** Literal coordinates, source identity, opposite strands and a physical end
	 * test the cache independently of any offset or network implementation. */
	private static void geometryRows(final File root, final NcrnaFamily family, final String member){
		final LineParser1 p=new LineParser1('\t');
		int count=0;
		boolean complete=false;
		try(Input input=new Input(path(root, "result.geometry.tsv"))){
			check(new String(input.reader.nextLine(), java.nio.charset.StandardCharsets.UTF_8).equals(EvaluateNcrnaRealEndpoints.GEOMETRY_HEADER),
				"Geometry schema must preserve explicit source/model and physical bounds");
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				if(row[0]=='#'){
					check(p.terms()==2 && p.termEquals("#COMPLETE", 0) && p.parseInt(1)==7 && count==7
						&& input.reader.nextLine()==null, "Exactly seven matched calls and one terminal marker are required");
					complete=true;
					break;
				}
				final int i=count++;
				final int model=i==6 ? 1 : 0;
				final int length=family.library[model].length;
				final int start=i==2 ? 0 : i==4 ? 2 : 32;
				final int stop=start+length-1;
				check(p.terms()==16 && p.parseInt(0)==i+1 && p.termEquals(i<5 ? "A" : "B", 1)
					&& p.termEquals(i<5 ? "assemblyA" : "assemblyB", 2) && p.termEquals(member, 4)
					&& p.termEquals(family.modelNames[model], 5) && p.parseInt(8)==length+64
					&& p.parseInt(9)==(i==2 ? 1 : 0), "Source-qualified identity and exact GFF fixed bits must survive geometry extraction");
				check(p.parseInt(10)==start && p.parseInt(11)==stop
					&& p.parseInt(12)==(i==2 || i==4 || i==6 ? start : start-5)
					&& p.parseInt(13)==(i==4 || i==6 ? stop : stop-3)
					&& p.parseInt(14)==family.minLen && p.length(15)==0 && family.maxLen==Integer.MAX_VALUE,
					"Both strands must retain literal gene-sense raw/truth coordinates and the selected family's length bounds");
			}
		}
		check(complete && count==7 && !new File(path(root, "result.features.tsv")).exists()
			&& !new File(path(root, "result.scores.tsv")).exists(), "Geometry-only mode neither drops calls nor materializes neural inference");
	}
	private static int maximum(final String text, final boolean legacy){
		final LineParser1 parser=new LineParser1('\t');
		parser.set(("left\t"+text+"\tright").getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		return NcrnaTextIO.maximumLength(parser, 1, legacy);
	}

	/** Actual adapter path, using a closed geometry cache and the supplied frozen
	 * first-stage pair. Whole assemblies select5/2/0rows; physical and missing
	 * holds survive. Literal audit identities independently check its labels and
	 * center-frame appendages, including the strand/source and footer contracts. */
	private static void gateFixture(final File root, final String[] original, final String familyName,
			final String offsets, final String roster, final String list){
		final String folds=path(root, "gate_folds.tsv");
		write(folds, new ByteBuilder(prok.NcrnaEndpointFoldInputs.ASSEMBLY_HEADER).nl()
			.append("assemblyA\t").append(path(root, "a.fa")).append("\tINTERNAL_FOLD\t0\n")
			.append("assemblyB\t").append(path(root, "b.fa")).append("\tINTERNAL_FOLD\t1\n"));
		final ArrayList<String> args=new ArrayList<String>();
		for(final String key:new String[]{"familyroot", "ref", "hbm", "seeds", "newtables", "newnets", "seedlabels"}){
			args.add(key+"="+value(original, key));
		}
		args.addAll(Arrays.asList("mode=gatevectors", "family="+familyName, "in="+path(root, "result.geometry.tsv"),
			"folds="+folds, "sources="+roster, "sourceformat=native", "fastalist="+list, "modeloffsets="+offsets,
			"inputs=2", "expected=7", "panel=finite", "t=1", "unpigz=f", "pigz=f"));
		for(int fold=0; fold<=2; fold++){
			final String prefix=path(root, "gate"+fold);
			args.add("fold="+fold);
			args.add("out="+prefix);
			new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process();
			for(final String seed:value(original, "seedlabels").split(",")){
				int proposed=0;
				for(final String end:new String[]{"5prime", "3prime"}){
					proposed+=gateAudit(prefix+".s"+seed+"."+end, fold, end, fold==0 ? 5 : fold==1 ? 2 : 0);
				}
				check(fold==2 ? proposed==0 : proposed>0,
					"Each nonempty supplied fold/seed must exercise actual35-input proposal rows; an all-raw-hold adapter cannot pass vacuously");
			}
		}
		args.add("fold=-1");
		args.add("out="+path(root, "gate_bad_outer"));
		fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process());
		check(!new File(path(root, "gate_bad_outer.features.tsv")).exists(), "Outer mode cannot silently discard internally assigned calls");
		final String outerFolds=path(root, "gate_outer_folds.tsv");
		write(outerFolds, new ByteBuilder(prok.NcrnaEndpointFoldInputs.ASSEMBLY_HEADER).nl()
			.append("assemblyA\t").append(path(root, "a.fa")).append("\tOUTER_TEST\t-1\n")
			.append("assemblyB\t").append(path(root, "b.fa")).append("\tOUTER_TEST\t-1\n"));
		args.add("folds="+outerFolds);
		args.add("out="+path(root, "gate_outer"));
		new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process();
		for(final String seed:value(original, "seedlabels").split(",")){
			int proposed=0;
			for(final String end:new String[]{"5prime", "3prime"}){proposed+=gateAudit(path(root, "gate_outer.s"+seed+"."+end), -1, end, 7);}
			check(proposed>0, "The complete outer fixture must also exercise real proposal rows");
		}
		final ByteBuilder truncated=new ByteBuilder();
		try(Input input=new Input(path(root, "result.geometry.tsv"))){
			for(byte[] row; (row=input.reader.nextLine())!=null;){if(row[0]!='#'){truncated.append(row).nl();}}
		}
		write(path(root, "gate_truncated.geometry.tsv"), truncated);
		args.add("in="+path(root, "gate_truncated.geometry.tsv"));
		args.add("out="+path(root, "gate_bad_truncated"));
		fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process());
		check(!new File(path(root, "gate_bad_truncated.features.tsv")).exists(), "A complete-looking cache without its terminal marker cannot become gate rows");
		geometryMaximum(root, "gate_legacy.geometry.tsv", true, "2147483647");
		args.add("in="+path(root, "gate_legacy.geometry.tsv"));
		args.add("out="+path(root, "gate_legacy"));
		new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process();
		for(final String seed:value(original, "seedlabels").split(",")){
			for(final String end:new String[]{"5prime", "3prime"}){
				equalRows(path(root, "gate_outer.s"+seed+"."+end+".vectors.tsv.gz"), path(root, "gate_legacy.s"+seed+"."+end+".vectors.tsv.gz"));
				gateAudit(path(root, "gate_legacy.s"+seed+"."+end), -1, end, 7);
			}
		}
		geometryMaximum(root, "gate_bad_max.geometry.tsv", false, "4294967496");
		args.add("in="+path(root, "gate_bad_max.geometry.tsv"));
		args.add("out="+path(root, "gate_bad_max"));
		fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])).process());
		check(!new File(path(root, "gate_bad_max.features.tsv")).exists(), "Overflowed maximum text must fail before gate output");
		args.add("fold=0");
		args.add("gff="+path(root, "must_not_reread.gff"));
		fails(()->new EvaluateNcrnaRealEndpoints(args.toArray(new String[args.size()])));
	}

	/** Literal old and malformed cache controls, independent of the maximum parser. */
	private static void geometryMaximum(final File root, final String name, final boolean legacy, final String maximum){
		final LineParser1 parser=new LineParser1('\t');
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path(root, "result.geometry.tsv")); Output output=new Output(path(root, name))){
			input.reader.nextLine();
			output.writer.println(legacy ? EvaluateNcrnaRealEndpoints.LEGACY_GEOMETRY_HEADER : EvaluateNcrnaRealEndpoints.GEOMETRY_HEADER);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				text.clear();
				if(line[0]=='#'){text.append(line);}
				else{
					parser.set(line);
					parser.setBounds(15);
					text.append(line, 0, parser.a()).append(maximum);
				}
				output.writer.print(text.nl());
			}
		}
	}

	private static int gateAudit(final String prefix, final int fold, final String end, final int expected){
		final LineParser1 p=new LineParser1('\t'), v=new LineParser1('\t'), s=new LineParser1('\t');
		int count=0, proposals=0, inferred=0;
		boolean complete=false;
		try(Input audit=new Input(prefix+".rows.tsv.gz"); Input matrix=new Input(prefix+".vectors.tsv.gz");
				Input scores=new Input(prefix+".scores.tsv.gz")){
			check(Arrays.equals(matrix.reader.nextLine(), "#dims\t35\t1\t0".getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"Real gate training matrix has its separately versioned35-input shape");
			byte[] scoreLine=scores.reader.nextLine();
			while(scoreLine!=null && scoreLine[0]=='#'){scoreLine=scores.reader.nextLine();}
			check(scoreLine!=null && new String(scoreLine, java.nio.charset.StandardCharsets.UTF_8).startsWith("row1\tstatus\t"),
				"Saved first-stage scores need their explicit raw-bit header");
			for(byte[] row; (row=audit.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("#COMPLETE", 0)){
					check(p.terms()==3 && p.parseInt(1)==expected && count==expected && p.parseInt(2)==proposals
						&& audit.reader.nextLine()==null && matrix.reader.nextLine()==null, "Every selected call and proposed gate row must close exactly once");
					s.set(scores.reader.nextLine());
					check(s.terms()==3 && s.termEquals("#COMPLETE", 0) && s.parseInt(1)==expected && s.parseInt(2)==inferred
						&& scores.reader.nextLine()==null, "Scores close at all selected calls and actual single-pass inference count");
					complete=true;
					break;
				}
				if(row[0]=='#' || p.termEquals("row1", 0)){continue;}
				count++;
				check(p.terms()==28 && (fold==-1 || p.termEquals(fold==0 ? "assemblyA" : "assemblyB", 2)) && p.termEquals(end, 20),
					"Gate rows retain every field and only the explicitly held assembly");
				final int id=p.parseInt(0);
				check(fold==-1 ? id>=1 && id<=7 : fold==0 ? id>=1 && id<=5 : id>=6 && id<=7,
					"Gate identity is the original cache row, not a synthetic trial");
				final int e=end.equals("5prime") ? 0 : 1;
				final int raw=p.parseInt(10+e), truth=p.parseInt(12+e), center=p.parseInt(16+e);
				final long proposal=p.parseLong(23);
				s.set(scores.reader.nextLine());
				check(s.parseInt(0)==id && s.termEquals(p.parseString(21), 1), "Audit and first-stage score identities/statuses must agree");
				if(s.terms()==24){inferred++;}else{check(s.terms()==2, "Held ends do not fabricate first-stage scores");}
				check(p.parseLong(25)==(long)center-raw, "Applied offset retains actual gene-sense displacement");
				if(p.parseInt(22)>0){
					check(p.parseInt(22)==++proposals && proposal!=raw && p.termEquals("PROPOSAL", 21) && s.terms()==24,
						"Only distinct-from-raw scored candidates enter the dense gate matrix");
					v.set(matrix.reader.nextLine());
					final boolean better=Math.abs(proposal-truth)<Math.abs((long)raw-truth);
					check(v.terms()==36 && v.parseInt(35)==(better ? 1 : 0) && p.parseInt(26)==(better ? 1 : 0),
						"Strict raw-to-proposal labels must agree independently in matrix and audit");
					for(int i=0; i<22; i++){
						check(Float.floatToRawIntBits(Float.parseFloat(v.parseString(i)))==s.parseInt(2+i), "Gate's first22inputs preserve the saved first-stage evidence");
					}
					check(v.parseFloat(32)==Math.abs((long)center-raw) && v.parseFloat(33)==(long)center-raw
						&& v.parseFloat(34)==(p.parseInt(24)==0 ? 0 : 1), "The three centered appendages have distinct exact meanings");
				}else{check(proposal==raw && p.parseInt(26)==0, "Every deterministic hold rejects to RAW, never to the median center");}
				if(id==3 && e==0){check(p.termEquals("FIXED", 21) && proposal==raw, "Physical fixed starts never enter gate inference");}
				if(id==4){check(p.termEquals("MISSING_FEATURE_HOLD", 21) && proposal==raw && s.terms()==2,
					"An ambiguous word at each centered endpoint must hold RAW without invented features or first-stage scores");}
				if(id==5){check(p.termEquals("INVALID_CENTER_RAW_HOLD", 21) && proposal==raw, "Invalid center preserves raw without inference");}
			}
		}
		check(complete, "An interrupted gate audit cannot pass as complete");
		return proposals;
	}

	private static String option(final String[] args, final String name, final String fallback){
		String result=fallback;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair[0].equalsIgnoreCase(name)){result=pair.length>1 ? pair[1] : null;}
		}
		return result;
	}
	private static String value(final String[] args, final String name){
		final String result=option(args, name, null);
		check(result!=null, "Missing finite resource binding: "+name);
		return result;
	}
	private static void write(final String path, final ByteBuilder text){try(Output out=new Output(path)){out.writer.print(text);}}
	private static void equalRows(final String first, final String second){
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] row=a.reader.nextLine(); row!=null; row=a.reader.nextLine()){
				check(Arrays.equals(row, b.reader.nextLine()), "Native and reordered legacy rosters must preserve identical source-qualified assay outputs");
			}
			check(b.reader.nextLine()==null, "Legacy adapter output must close at the same row count");
		}
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void fails(final Runnable action){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Malformed provenance must fail loudly");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
