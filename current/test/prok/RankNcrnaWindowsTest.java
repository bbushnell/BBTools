package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.CallGenes;
import prok.CallGenesIteration;
import prok.ConservedRnaSeedIndex;
import prok.GeneCaller;
import prok.GeneModel;
import prok.GradeNcrnaReference;
import prok.NcrnaDdlOracle;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaMemberLengthBounds;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import prok.RankNcrnaWindows;
import shared.Shared;
import shared.Timer;
import stream.Read;
import structures.ByteBuilder;

/** Full native-caller observation parity and complete report-denominator fixture.
 * Synthetic genomes embed one real consensus sequence, while synthetic one-member
 * profiles preserve the complete native model order. These are not construction
 * resources or a performance/accuracy panel. The control observes ordinary attempts
 * without doing oracle/ranking work; final GFF and every attempt must match bytes.
 * @author Brian Bushnell, Raiden
 */
public final class RankNcrnaWindowsTest{

	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, RankNcrnaWindowsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown full-ranking fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply out=fresh_directory without inputs, overwrite or sampling");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "The native parity fixture owns a fresh output directory");
			new RankNcrnaWindowsTest(root).run();
			pp.outstream.println("NCRNA_DDL_RANK_DRIVER_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private RankNcrnaWindowsTest(final File root_){root=root_;}

	private void run() throws IOException{
		checkQuantiles();
		final NcrnaFamily family=NcrnaFamilyLoader.load("cobalamin_rf00174");
		final ConservedRnaSeedIndex seeds=new ConservedRnaSeedIndex(new map.LongHashSet[]{family.kmerSet});
		byte[] gene=null;
		for(final byte[] reference:family.library){
			if(reference.length>=family.minLen && seeds.scan(reference).hits(0).length>0){
				gene=reference;
				break;
			}
		}
		check(gene!=null, "The real library needs one seed-supported consensus for this native fixture");
		final byte[] genome=new byte[gene.length+800];
		Arrays.fill(genome, (byte)'N');
		System.arraycopy(gene, 0, genome, 200, gene.length);
		write("genome.fa", new ByteBuilder(">unique description\n").append(genome)
			.append("\n>ambiguous description\n").append(genome).nl());
		write("zero.fa", new ByteBuilder(">zero description\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"));
		metadata(gene.length);
		members(family);
		final String config="in="+path("genome.fa")+','+path("zero.fa")+
			"\ncobalamin=t\ntax=f\npasses=1\nordered=t\nt=1\n";
		write("caller.config", new ByteBuilder(config));
		final ByteBuilder controlAttempts=new ByteBuilder();
		final long[] observedWindows={0};
		final CallGenes control=new CallGenes(Parser.parseConfig(new String[]{"config="+path("caller.config"),
			"out="+path("control.gff"), "outstats=null", "ow=f"}));
		control.run(new Timer(), ()->control.setIteration((input, caller, model)->{
			final int slot=input.equals(path("genome.fa")) ? 0 : input.equals(path("zero.fa")) ? 1 : -1;
			check(slot>=0, "The ordinary control must see only the two literal genomes");
			return new Control(slot, caller, model, controlAttempts, observedWindows);
		}));
		RankNcrnaWindows.main(arguments("caller.config", "rank", family.library.length));
		check(observedWindows[0]>0 && controlAttempts.length()>0, "The whole-caller control must genuinely capture and attempt models");
		check(Arrays.equals(Files.readAllBytes(new File(path("control.gff")).toPath()),
			Files.readAllBytes(new File(path("rank.calls.gff")).toPath())),
			"The synchronous oracle/profile reports must preserve the entire ordinary gene-path GFF");
		final ByteBuilder actualAttempts=new ByteBuilder();
		try(Input in=new Input(path("rank.attempts.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){actualAttempts.append(row).nl();}
		}
		check(Arrays.equals(controlAttempts.toBytes(), actualAttempts.toBytes()),
			"Every actual native attempt must preserve event, input, strand-window binding, rank, score and initial-pass flag");
		checkReports(observedWindows[0], family.library.length);
		checkNativeOrder(family.library.length);
		checkCallerSelection(config, family.library.length);
	}

	private String[] arguments(final String caller, final String out, final int models){
		return new String[]{"caller="+path(caller), "out="+path(out),
			"reference="+path("reference.tsv"), "familymap="+path("families.tsv"),
			"subtypes="+path("subtypes.tsv"), "witnesses="+path("witnesses.tsv"),
			"events="+path("events.tsv"), "members="+path("members.fa"), "ledger="+path("ledger.tsv"),
			"assignments="+path("assignments.tsv"), "census="+path("models.tsv"),
			"membercounts="+path("membercounts.tsv"), "expectedmodels="+models,
			"ks=5", "widths=32", "t=1"};
	}

	/** Exercise the expanded caller file, not just the ranking command's own flags. */
	private void checkCallerSelection(final String config, final int models){
		final String[] flags={"reads=1", "maxreads=1", "samplerate=0.5"};
		for(int i=0; i<flags.length; i++){
			final String name="limited_"+i, output=name+"_result";
			write(name+".config", new ByteBuilder(config).append(flags[i]).nl());
			boolean rejected=false;
			try{RankNcrnaWindows.main(arguments(name+".config", output, models));}
			catch(IllegalArgumentException expected){
				rejected=expected.getMessage().startsWith("Caller configuration must not limit or sample contigs");
			}
			check(rejected, "A read limit or sample rate hidden inside caller.config must fail before outputs: "+flags[i]);
			final File[] files=root.listFiles();
			check(files!=null, "The owned fixture directory must remain readable");
			for(final File file:files){
				check(!file.getName().startsWith(output+'.'), "Invalid caller selection must not create a partial ranking report");
			}
		}
	}

	/** Explicit native control, with no exhaustive verifier, sketches or truth join. */
	private final class Control implements CallGenesIteration.Worker{
		Control(final int slot_, final GeneCaller caller_, final GeneModel model_,
				final ByteBuilder attempts, final long[] windows){
			slot=slot_;
			caller=caller_;
			model=model_;
			int target=-1;
			for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
				if(GeneCaller.ncrnaFamilies.get(i).name.equals("cobalamin_rf00174")){target=i;}
			}
			check(target>=0, "The whole-caller control must use the native RF00174 worker");
			caller.ncrnaScavenger(target).setWindowCaptureSink((name, strand, pass, start, stop, source, bases, window)->{
				event++;
				windows[0]++;
			});
			caller.ncrnaScavenger(target).setModelAttemptSink((name, strand, pass, start, stop, rank, m, shared, identity, length, initial)->{
				attempts.append(slot).tab().append(slot==0 ? "I" : "ZERO").tab().append(event).tab()
					.append(contig).tab().append(NcrnaDdlOracle.seqid(name)).tab().append(rank).tab()
					.append(m).tab().append(shared).tab().append(Float.floatToRawIntBits(identity)).tab()
					.append(length).tab().append(initial).nl();
			});
		}
		@Override
		public ArrayList<Orf> call(final Read read){
			contig++;
			final byte[] before=read.bases.clone();
			final ArrayList<Orf> result=caller.callGenes(read, model, true);
			check(Arrays.equals(before, read.bases), "The complete ordinary control restores each original contig");
			return result;
		}
		@Override
		public void finish(){}
		private final int slot;
		private final GeneCaller caller;
		private final GeneModel model;
		private long event=-1, contig=-1;
	}

	/** Five physical parents: three targets, two competitors; empty evidence uses explicit priors. */
	private void metadata(final int length){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> parents=new ArrayList<Row>(), zero=new ArrayList<Row>();
		parents.add(codec.create("I", new Locus("unique", "cobalamin", 201, 200+length, '+'),
			false, null, "RF00174", path("witnesses.tsv")));
		uniqueParent=parents.get(0).id;
		parents.add(codec.create("I", new Locus("ambiguous", "cobalamin", 201, 200+length, '+'),
			false, null, "RF00174", path("witnesses.tsv")));
		ambiguousParent=parents.get(1).id;
		parents.add(codec.create("I", new Locus("ambiguous", "cobalamin", 202, 201+length, '+'),
			false, null, "RF01482", path("witnesses.tsv")));
		otherRight=parents.get(2).id;
		// Two opposing one-base shifts ensure every span compatible with the target
		// also overlaps a scored competitor, without duplicate physical parent IDs.
		parents.add(codec.create("I", new Locus("ambiguous", "cobalamin", 200, 199+length, '+'),
			false, null, "RF01482", path("witnesses.tsv")));
		otherLeft=parents.get(3).id;
		zero.add(codec.create("ZERO", new Locus("zero", "cobalamin", 3, 25, '+'),
			false, null, "RF00174", path("witnesses.tsv")));
		codec.write(path("reference.gff"), "I", parents);
		codec.write(path("zero.gff"), "ZERO", zero);
		write("reference.tsv", new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl()
			.append("I\t").append(path("genome.fa")).tab().append(path("reference.gff")).nl()
			.append("ZERO\t").append(path("zero.fa")).tab().append(path("zero.gff")).nl());
		write("families.tsv", new ByteBuilder(GradeNcrnaReference.RULE_HEADER).nl()
			.append("cobalamin\tcobalamin_RF00174\tRF00174\ttrue\ncobalamin\tcobalamin_RF01482\tRF01482\ttrue\n"));
		write("subtypes.tsv", new ByteBuilder("subtype\tmodel_prefix\treference_model\n")
			.append("main\tcobalamin_RF00174\tRF00174\nother\tcobalamin_RF01482\tRF01482\n"));
		write("events.tsv", new ByteBuilder(EVENT_HEADER).nl());
		final ByteBuilder witnesses=new ByteBuilder("annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation\n");
		parents.addAll(zero);
		for(final Row row:parents){
			BuildNcrnaReference.field(witnesses, row.id).tab().append(row.input)
				.append("\tPROTECTED_PRIOR\t.\t.\t.\t").append(row.model).tab()
				.append(row.locus.start).tab().append(row.locus.stop).append("\t2\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n");
		}
		write("witnesses.tsv", witnesses);
	}

	/** One synthetic kept member per native model; no public construction census is relabelled. */
	private void members(final NcrnaFamily family){
		final ByteBuilder fasta=new ByteBuilder();
		final ByteBuilder ledger=new ByteBuilder("locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header\n");
		final ByteBuilder assignments=new ByteBuilder("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n");
		final ByteBuilder census=new ByteBuilder(NcrnaMemberLengthBounds.HEADER).nl();
		for(int m=0; m<family.library.length; m++){
			final byte[] sequence=family.library[m].clone();
			for(int i=0; i<sequence.length; i++){
				if(AminoAcid.baseToNumber[sequence[i]]<0){sequence[i]='A';}
			}
			final String id="RF00174_fixture_"+m;
			final String raw=id+" source=tid_1_species status=ANNOTATION_FULL_COPY";
			final String header=raw+" member_family=RF00174 build_group=RF00174";
			final String name=NcrnaDdlOracle.seqid(family.modelNames[m]);
			fasta.append('>').append(header).nl().append(sequence).nl();
			ledger.append(id).append("\tRF00174\tRF00174\tRF00174\tACCEPTED_RF00174_PRODUCT\tNA\t")
				.append(sequence.length).tab().append(raw).nl();
			assignments.append(header).append("\t1\tunknown\tRF00174_fixture\t").append(name).append("\tkept\n");
			census.append(name).append("\tRF00174\t").append(sequence.length).append("\t1");
			for(int i=0; i<5; i++){census.tab().append(sequence.length);}
			census.append("\t1\tPURE\t1\t1\t0\t0\t0\t0\t0\t0\n");
		}
		write("members.fa", fasta);
		write("ledger.tsv", ledger);
		write("assignments.tsv", assignments);
		write("models.tsv", census);
		final int n=family.library.length;
		write("membercounts.tsv", new ByteBuilder("emitted_family\tinput_members\toriginal_kept\toriginal_small\toriginal_orphan\tadded_memberships\tconstruction_memberships\tmodels\tpure\tmixed\tempty\tpurity_cut\tgroup_min\tgroup_p10\tgroup_median\tgroup_p90\tgroup_max\n")
			.append("RF00174\t").append(n).tab().append(n).append("\t0\t0\t0\t").append(n).tab()
			.append(n).tab().append(n).append("\t0\t0\t0.8\t17\t17\t17\t17\t17\n"));
	}

	private void checkReports(final long windows, final int models){
		final LineParser1 p=new LineParser1('\t');
		int inputs=0, parents=0, loci=0, zeroLoci=0, ambiguousLoci=0, summary=0, geometryRows=0;
		long rankingCpu=0, geometryCpu=0, profileCpu=-1;
		try(Input in=new Input(path("rank.inputs.tsv"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==13 && p.parseInt(7)==(inputs==0 ? 2 : 1), "Input target denominators retain the zero-window genome");
				check(p.parseLong(5)==(inputs==0 ? windows : 0), "Input window totals match the independent ordinary observer");
				rankingCpu+=p.parseLong(11);
				inputs++;
			}
		}
		try(Input in=new Input(path("rank.parents.tsv.gz"))){
			in.reader.nextLine();
			while(in.reader.nextLine()!=null){parents++;}
		}
		try(Input in=new Input(path("rank.loci.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==12, "Every target/method row retains complete missing-rank and reachability fields");
				if(p.termEquals("ZERO", 1)){
					check(p.parseLong(4)==0 && p.parseInt(8)==-1 && p.termEquals("NO_CAPTURED_WINDOW", 9),
						"A target without a seed/window remains missing in every method, never rank one");
					zeroLoci++;
				}
				final int rank=p.parseInt(8);
				check(rank==-1 || rank>0 && rank<=models, "Every reported rank lies in its actual native model population");
				if(p.termEquals(uniqueParent, 2) && p.parseInt(3)>0){
					check(rank>0, "Each DDL method must rank an acceptable model for the isolated exact-consensus parent");
				}
				if(p.termEquals(ambiguousParent, 2)){
					check(p.parseLong(5)>0 && p.parseLong(6)==p.parseLong(5) && p.parseLong(7)==0
						&& rank==-1 && p.termEquals("AMBIGUOUS_SPANS_ONLY", 9)
						&& p.termEquals("false", 10) && p.termEquals("false", 11),
						"Other-subtype scored parents must prevent every ambiguous target span from receiving a unique rank or reachability credit");
					ambiguousLoci++;
				}
				loci++;
			}
		}
		try(Input in=new Input(path("rank.summary.tsv"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==24 && p.parseInt(4)==3 && p.parseInt(7)>=1,
					"Every distribution/cutoff row retains all three targets and the explicit no-window locus");
				check(p.parseLong(11)+p.parseLong(12)+p.parseLong(13)==p.parseLong(6),
					"Improved/tied/worsened comparisons conserve the common native-rankable denominator");
				summary++;
			}
		}
		try(Input in=new Input(path("rank.geometry_cpu.tsv"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==8 && p.parseInt(1)==5 && p.parseInt(2)==32 && p.parseLong(3)>=0
					&& p.parseLong(4)>=0 && p.parseLong(5)>=0 && p.parseLong(6)==(long)models*32*2
					&& p.parseLong(7)==(long)models*32*4,
					"Each requested geometry retains its own measured CPU scopes and explicitly primitive profile payload sizes");
				geometryRows++;
				profileCpu=p.parseLong(3);
				geometryCpu+=p.parseLong(4)+p.parseLong(5);
			}
		}
		boolean profileClock=false, memberClock=false;
		try(Input in=new Input(path("rank.work.tsv"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(p.termEquals("profile_cpu_ns", 0)){profileClock=p.parseLong(1)==profileCpu;}
				if(p.termEquals("member_ingestion_owner_cpu_ns", 0)){memberClock=p.parseLong(1)>=0;}
			}
		}
		check(profileClock && memberClock && rankingCpu>=geometryCpu,
			"The one-geometry profile clock must match its total, member loading stays separate, and geometry costs fit whole-worker ranking CPU");
		check(inputs==2 && parents==5 && loci==12 && zeroLoci==4 && ambiguousLoci==4 && summary==88 && geometryRows==1,
			"Two inputs/five parents/three targets/four methods/two distributions/eleven cutoffs and one geometry must conserve exactly");
	}

	/** Reconstructs the baseline order from complete emitted counts and independent
	 * compatible-parent links, rather than invoking the production ranking helper. */
	private void checkNativeOrder(final int models){
		final ObjectMap<String, NativeWindow> windows=new ObjectMap<String, NativeWindow>(String.class, NativeWindow.class);
		final ArrayList<NativeWindow> observed=new ArrayList<NativeWindow>();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path("rank.oracle.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				final String key=p.parseInt(0)+"_"+p.parseLong(2);
				NativeWindow window=windows.get(key);
				if(window==null){
					window=new NativeWindow(models);
					windows.put(key, window);
					observed.add(window);
				}
				final int model=p.parseInt(5);
				check(model>=0 && model<models && !window.seen[model], "Every event/model must occur once in the complete oracle ledger");
				window.seen[model]=true;
				window.counts[model]=p.parseInt(7);
				window.rows++;
			}
		}
		try(Input in=new Input(path("rank.links.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				final NativeWindow window=windows.get(p.parseInt(0)+"_"+p.parseLong(2));
				check(window!=null, "A parent link must bind an existing exhaustive event");
				final int model=p.parseInt(5);
				if(p.termEquals(uniqueParent, 6) && p.termEquals("true", 9)){window.acceptable[model]=true;}
				if(p.termEquals(ambiguousParent, 6)){
					check(p.parseInt(8)>1 && p.termEquals("false", 9), "Each ambiguous target link must retain multiple scored parents without unique attribution");
					window.ambiguous[model]=true;
				}
				if(p.termEquals(otherLeft, 6) || p.termEquals(otherRight, 6)){window.other[model]=true;}
			}
		}
		int checked=0;
		try(Input in=new Input(path("rank.ranks.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(p.parseInt(6)!=0 || !p.termEquals(uniqueParent, 5)){continue;}
				final NativeWindow window=windows.get(p.parseInt(0)+"_"+p.parseLong(2));
				check(window!=null && window.rows==models, "A native rank must retain its entire original model population");
				int selected=-1;
				for(int m=0; m<models; m++){
					if(window.acceptable[m] && window.counts[m]>0
						&& (selected<0 || window.counts[m]>window.counts[selected])){selected=m;}
				}
				int rank=-1, first=-1, last=-1;
				if(selected>=0){
					rank=first=1;
					last=0;
					for(int m=0; m<models; m++){
						if(window.counts[m]>window.counts[selected]){
							rank++;
							first++;
							last++;
						}else if(window.counts[m]==window.counts[selected]){
							last++;
							if(m<selected){rank++;}
						}
					}
				}
				check(p.parseInt(11)==rank && p.parseInt(12)==first && p.parseInt(13)==last && p.parseInt(14)==selected,
					"The emitted baseline rank/tie must match descending native counts and original-index tie order");
				checked++;
			}
		}
		int ambiguousLinks=0;
		for(final NativeWindow window:observed){
			check(window.rows==models, "Every captured event must retain every native model, including rejected and zero-support models");
			for(int model=0; model<models; model++){
				if(!window.ambiguous[model]){continue;}
				check(window.other[model], "An ambiguous target model must retain its compatibility link to a scored competing-subtype parent");
				ambiguousLinks++;
			}
		}
		check(ambiguousLinks>0, "The whole-driver ambiguity fixture must actually produce multi-parent compatibility links");
		check(checked>0, "The native-order cross-check must exercise the isolated physical parent");
	}
	private static final class NativeWindow{
		NativeWindow(final int models){
			counts=new int[models];
			seen=new boolean[models];
			acceptable=new boolean[models];
			ambiguous=new boolean[models];
			other=new boolean[models];
		}
		final int[] counts;
		final boolean[] seen, acceptable, ambiguous, other;
		int rows;
	}

	private static void checkQuantiles(){
		final long[] histogram={0, 2, 0, 1, 0, 1};
		check(RankNcrnaWindows.quantile(histogram, .5)==1
			&& RankNcrnaWindows.quantile(histogram, .9)==5 && RankNcrnaWindows.quantile(histogram, 1)==5,
			"The literal ranks1,1,3,5 have nearest-rank median1 and upper quantiles5");
		check(RankNcrnaWindows.quantile(new long[4], .5)==-1, "An empty rank denominator has no numeric quantile");
	}
	private void write(final String file, final ByteBuilder bytes){try(Output out=new Output(path(file))){out.writer.print(bytes);}}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private String uniqueParent, ambiguousParent, otherLeft, otherRight;
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
}
