package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.IntObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.LongList;

/** Reduces disjoint completed saved-window shards from their per-gene winners.
 * Uses A_SampleBasic2's native parser/output lifecycle and byte metadata inputs.
 * Every class2 gene must have exactly one owning shard and every stage/condition
 * exactly once. Curves and histograms are recomputed from those winners rather
 * than adding partial statistics or counting a gene once per overlapping window.
 * No sequence, model or alignment is read or executed here.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaRejectedWindows{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaRejectedWindows.class, false);
		try{
			new SummarizeNcrnaRejectedWindows(pp.args).process();
			pp.outstream.println("NCRNA_REJECTED_WINDOW_SUMMARY_COMPLETE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeNcrnaRejectedWindows(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("misses")){misses=b;}
			else if(a.equals("population")){population=b;}
			else if(a.equals("expectedlibraries")){expectedLibraries=Parse.parseIntKMG(b);}
			else if(a.equals("expectedcaptures")){expectedCaptures=Parse.parseIntKMG(b);}
			else if(a.equals("expectedrows")){expectedRows=Parse.parseIntKMG(b);}
			else if(a.equals("expectedoriginals")){expectedOriginals=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown attempt-summary argument: "+arg);}
		}
		manifest=parser.in1;
		prefix=parser.out1;
		require(manifest!=null && prefix!=null && misses!=null && population!=null && expectedLibraries>0
			&& expectedCaptures>0 && expectedRows>0 && expectedOriginals>0 && !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=completed_shards.tsv misses=full_gene_misses population=prefix out=fresh and explicit population counts without sampling");
		RrnaResourceIO.fresh(prefix, ".curves.tsv", ".histograms.tsv", ".own_models.tsv", ".work.tsv");
		Shared.setThreads(1);
	}
	private void process(){
		loadManifest();
		loadExpected();
		loadOriginals();
		for(final Capture capture:captures){readCapture(capture);}
		for(final Library library:libraries){
			for(final Gene gene:library.order){require(gene.mask==ALL_CONDITIONS, "Every class2 original must occur in exactly one shard with every stage/condition; gene="+gene.id);}
		}
		report();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Population Joins       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadManifest(){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(Input in=new Input(manifest)){
			RrnaResourceIO.header(in.reader, "library\tcapture\tresult_prefix", manifest);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==3, "Every completed shard requires its library, namespace and result prefix");
				final String name=p.parseString(0), namespace=p.parseString(1), source=p.parseString(2);
				Library library=null;
				for(final Library candidate:libraries){if(candidate.name.equals(name)){library=candidate;break;}}
				if(library==null){library=new Library(name);libraries.add(library);}
				require(seen.add(name+"/"+namespace) && source.startsWith("/"), "Captured namespaces must be unique within a library with absolute prefixes");
				captures.add(new Capture(library, namespace, source));
			}
		}
		require(libraries.size()==expectedLibraries && captures.size()==expectedCaptures, "The complete declared library/shard population must be present");
	}
	private void loadExpected(){
		int count=0;
		try(Input in=new Input(misses)){
			require(in.reader.nextLine()!=null, "The accepted all-original miss report needs its header");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("CAPTURED_NO_GENE_ACCEPTANCE", 4)){continue;}
				final Library library=library(0);
				require(library!=null, "A summary cannot omit a library represented in the accepted class2 population");
				final int id=p.parseInt(1);
				require(!library.genes.contains(id), "Expected class2 IDs must be unique per library");
				final Gene gene=new Gene(id, bool(3));
				library.genes.put(id, gene);
				library.order.add(gene);
				count++;
			}
		}
		require(count==expectedRows, "Class2 denominators must be the accepted full library-gene population, never the number reported by completed shards");
	}
	private void loadOriginals(){
		int count=0;
		try(Input in=new Input(population+".original_calls.tsv.gz")){
			RrnaResourceIO.header(in.reader, PrepareNcrnaRosterBaseline.ORIGINAL_HEADER+"\tbaseline_input_id\tbaseline_input_state", population);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				count++;
				for(final Library library:libraries){
					final Gene gene=library.genes.get(p.parseInt(0));
					if(gene==null){continue;}
					require(!gene.original && gene.incomplete==bool(21), "Class2 original coverage must have one unchanged native metadata binding");
					gene.original=true;
					gene.shortBoth=p.parseInt(18)-1==44 && p.parseInt(20)-p.parseInt(19)==20;
					gene.start=p.parseInt(4)-1;
					gene.stop=p.parseInt(5)-1;
				}
			}
		}
		require(count==expectedOriginals, "Original coordinate metadata cannot be sampled while aggregating class2 subsets");
		for(final Library library:libraries){
			for(final Gene gene:library.order){require(gene.original, "No class2 original may lose its coordinate/coverage binding");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Completed Winners      ----------------*/
	/*--------------------------------------------------------------*/

	private void readCapture(final Capture capture){
		int originals=0;
		long rows=0;
		try(Input in=new Input(capture.prefix+".best.tsv.gz")){
			RrnaResourceIO.header(in.reader, BEST_HEADER, capture.prefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==20 && p.termEquals(capture.library.name, 0) && p.termEquals(capture.namespace, 1), "A raw winner row cannot cross library or shard identities");
				final Gene gene=capture.library.genes.get(p.parseInt(2));
				require(gene!=null && gene.incomplete==bool(3) && gene.shortBoth==bool(4), "Reported original IDs and native subsets must match the independent complete class2 population");
				final int stage=term(5, STAGES), condition=term(6, CONDITIONS), bit=1<<(stage*10+condition);
				require((gene.mask&bit)==0 && (gene.capture==null || gene.capture==capture), "Every original has exactly one owning shard and one row per stage/condition");
				gene.capture=capture;
				gene.mask|=bit;
				if(stage==0 && condition==0){originals++;}
				final boolean known=!p.termEquals("NA", 7), own=!p.termEquals("NA", 17);
				if(gene.ownState<0){gene.ownState=own ? 1 : 0;}
				require(gene.ownState==(own ? 1 : 0), "An original gene's construction-model availability cannot change between stages or conditions");
				float identity=Float.NaN;
				int span=-1;
				boolean geometry=false;
				if(known){
					identity=p.parseFloat(7);
					span=p.parseInt(9);
					final int start=p.parseInt(10), stop=p.parseInt(11);
					geometry=qualifies(gene.start, gene.stop, start, stop);
					require(Float.isFinite(identity) && identity>=0 && identity<=1 && start>=0 && stop>=start && span==stop-start+1
						&& bool(13)==(identity<.64f) && bool(14)==(span<109) && bool(15)==!geometry,
						"Raw winner identity, inclusive span and gate flags must agree with original coordinates");
					final int winningStage=term(19, STAGES);
					require(winningStage<2 && (stage==2 || winningStage==stage), "A union winner must retain its actual primary/edge source stage");
					require(winningStage==1 ? bool(18)==(identity>=NcrnaScavenger.EDGE_IDENTITY_FLOOR) : p.termEquals("NA",18),
						"Only physical-edge attempts have the separately reported native edge floor");
					require(condition==1 || Math.min(gene.stop,stop)>=Math.max(gene.start,start), "Positive-overlap winners cannot be invented from an unrelated interval");
					if(condition==1 || condition==2 || condition==9){require(sameColumns(8,17), "An own-model winner must name that gene's actual construction model");}
					if(condition>=3){
						final int cell=condition==9 ? 4 : condition-3;
						require(geometry && identity>=CUTS[cell/2] && (cell%2!=0 || span>=109)
							&& (winningStage==0 || bool(16) && bool(18)), "An eligible condition needs its own identity/span/geometry and genuine-edge gates");
					}
				}else{
					for(int field=7;field<=16;field++){require(p.termEquals("NA",field), "An unavailable attempt cannot retain invented identity, span or gate flags");}
					require(p.termEquals("NA",18) && p.termEquals("NA",19), "An unavailable winner cannot invent an edge floor or stage");
				}
				if(condition==1 || condition==2 || condition==9){require(own || !known, "A missing construction model cannot acquire an invented own-model winner");}
				capture.library.observe(gene, stage, condition, known, own, identity, span);
				rows++;
			}
		}
		require(rows==(long)originals*30, "Each shard must retain all thirty stage/condition rows per class2 original");
		final long[] totals=new long[SUMMARY_FIELDS.length];
		java.util.Arrays.fill(totals,-1);
		try(Input in=new Input(capture.prefix+".summary.tsv")){
			RrnaResourceIO.header(in.reader, "metric\tcount", capture.prefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==2, "Completed work summaries require named scalar counters");
				final long value=p.parseLong(1);
				require(value>=0, "Work counters cannot be negative");
				for(int field=0; field<SUMMARY_FIELDS.length; field++){
					if(!p.termEquals(SUMMARY_FIELDS[field],0)){continue;}
					require(totals[field]<0, "A completed summary cannot repeat a work/population counter: "+SUMMARY_FIELDS[field]);
					totals[field]=value;
				}
			}
		}
		for(final long value:totals){require(value>=0, "Every claimed work/population counter must occur exactly once before any count is added");}
		require(totals[0]==capture.library.order.size(), "The full-library class2 counter must match the independently loaded accepted population");
		require(totals[1]==originals, "Global unique-gene collection must agree with each separately counted shard summary");
		capture.library.windows=Math.addExact(capture.library.windows,totals[2]);
		capture.library.primary=Math.addExact(capture.library.primary,totals[3]);
		capture.library.edge=Math.addExact(capture.library.edge,totals[4]);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Output           ----------------*/
	/*--------------------------------------------------------------*/

	private void report(){
		try(Output curve=new Output(prefix+".curves.tsv"); Output hist=new Output(prefix+".histograms.tsv");
			Output own=new Output(prefix+".own_models.tsv"); Output work=new Output(prefix+".work.tsv")){
			curve.writer.println("library\tcohort\tstage\tclass2_genes\tidentity_cut\tspan_floor\tmeasured_eligible_genes_NOT_FINAL_RECALL");
			hist.writer.println("library\tcohort\tstage\tquantity\tlower_inclusive\tupper_exclusive\tgenes");
			own.writer.println("library\tcohort\tstage\tclass2_genes\tno_construction_model\town_model_eligible064_span109\town_model_not_eligible064_span109\town_model_best_identity_below064\town_model_no_positive_overlap");
			work.writer.println("library\tclass2_genes\tselected_saved_windows\tprimary_model_alignments\tphysical_edge_model_alignments");
			for(final Library library:libraries){
				work.writer.print(bytes.clear().append(library.name).tab().append(library.order.size()).tab().append(library.windows).tab().append(library.primary).tab().append(library.edge).nl());
				for(int c=0; c<3; c++){
					for(int stage=0; stage<3; stage++){
						for(int cell=0; cell<6; cell++){
							curve.writer.print(key(library, c, stage).tab().append(library.denominators[c]).tab().append(CUTS[cell/2],2).tab()
								.append(cell%2==0 ? 109 : 0).tab().append(library.curves[c][stage][cell]).nl());
						}
						final long[] counts=library.own[c][stage];
						require(counts[0]+counts[1]+counts[2]==library.denominators[c], "Missing, eligible and ineligible own-model states must conserve each full cohort");
						own.writer.print(key(library,c,stage).tab().append(library.denominators[c]).tab().append(counts[0]).tab().append(counts[1]).tab().append(counts[2]).tab().append(counts[3]).tab().append(counts[4]).nl());
						for(int quantity=0; quantity<4; quantity++){
							final LongList bins=library.hist[c][stage][quantity];
							long total=library.absent[c][stage][quantity];
							for(int bin=0; bin<bins.size; bin++){
								if(bins.get(bin)==0){continue;}
								total+=bins.get(bin);
								key(library,c,stage).tab().append(QUANTITIES[quantity]).tab();
								if(quantity==1){bytes.append(bin).tab().append(bin+1);}
								else{bytes.append(bin/50.0,2).tab().append(bin==50 ? "INCLUSIVE_1" : Double.toString((bin+1)/50.0));}
								hist.writer.print(bytes.tab().append(bins.get(bin)).nl());
							}
							hist.writer.print(key(library,c,stage).tab().append(QUANTITIES[quantity]).append("\tUNAVAILABLE\tNA\t").append(library.absent[c][stage][quantity]).nl());
							require(total==library.denominators[c], "Every recomputed full-population histogram must conserve its original cohort including unavailable attempts");
						}
					}
				}
			}
		}
	}
	private ByteBuilder key(Library library, int cohort, int stage){return bytes.clear().append(library.name).tab().append(COHORTS[cohort]).tab().append(STAGES[stage]);}
	private Library library(int column){
		for(final Library library:libraries){if(p.termEquals(library.name, column)){return library;}}
		return null;
	}
	private int term(int column, String[] values){
		for(int i=0; i<values.length; i++){if(p.termEquals(values[i], column)){return i;}}
		throw new IllegalArgumentException("Unknown stage/condition in completed winner table");
	}
	private boolean bool(int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Saved gate/subset flags must be explicit");
		return p.termEquals("true", column);
	}
	private boolean sameColumns(int first, int second){
		p.setBounds(first);
		final int from=p.a(), length=p.b()-p.a();
		p.setBounds(second);
		if(length!=p.b()-p.a()){return false;}
		for(int i=0; i<length; i++){if(p.line()[from+i]!=p.line()[p.a()+i]){return false;}}
		return true;
	}
	private static boolean qualifies(int a, int b, int c, int d){
		final long n=Math.max(0, Math.min(b,d)-Math.max(a,c)+1);
		return 2*n>=b-(long)a+1 && 2*n>=d-(long)c+1;
	}
	private static void require(boolean valid, String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Library{
		Library(String name_){
			name=name_;
			for(int c=0; c<3; c++){
				for(int stage=0; stage<3; stage++){
					for(int quantity=0; quantity<4; quantity++){hist[c][stage][quantity]=new LongList();}
				}
			}
		}
		void observe(Gene gene, int stage, int condition, boolean known, boolean hasOwn, float identity, int span){
			assert(gene.original && stage>=0 && stage<3 && condition>=0 && condition<10) : "Only validated original/condition rows may update the full-population histograms";
			for(int c=0; c<3; c++){
				if(c==1 && !gene.incomplete || c==2 && !gene.shortBoth){continue;}
				if(stage==0 && condition==0){denominators[c]++;}
				if(condition>=3 && condition<=8 && known){curves[c][stage][condition-3]++;}
				if(condition==9){own[c][stage][!hasOwn ? 0 : known ? 1 : 2]++;}
				if(condition==1 && hasOwn && known && identity<.64f){own[c][stage][3]++;}
				if(condition==2 && hasOwn && !known){own[c][stage][4]++;}
				if(condition>2){continue;}
				final int quantity=condition==0 ? 0 : condition+1;
				if(!known){
					absent[c][stage][quantity]++;
					if(condition==0){absent[c][stage][1]++;}
				}
				else{
					int bin=0;
					while(bin<50 && identity>=(float)((bin+1)/50.0)){bin++;}
					hist[c][stage][quantity].increment(bin);
					if(condition==0){hist[c][stage][1].increment(span);}
				}
			}
		}
		final String name;
		final IntObjectMap<Gene> genes=new IntObjectMap<Gene>(16,Gene.class);
		final ArrayList<Gene> order=new ArrayList<Gene>();
		final long[] denominators=new long[3];
		final long[][][] curves=new long[3][3][6], own=new long[3][3][5], absent=new long[3][3][4];
		final LongList[][][] hist=new LongList[3][3][4];
		long windows, primary, edge;
	}
	private static final class Gene{
		Gene(int id_, boolean incomplete_){
			id=id_;
			incomplete=incomplete_;
		}
		final int id;
		final boolean incomplete;
		boolean original, shortBoth;
		int start, stop, mask;
		int ownState=-1;
		Capture capture;
	}
	private static final class Capture{
		Capture(Library library_, String namespace_, String prefix_){
			library=library_;
			namespace=namespace_;
			prefix=prefix_;
		}
		final Library library;
		final String namespace, prefix;
	}
	private String manifest, misses, population, prefix;
	private int expectedLibraries, expectedCaptures, expectedRows, expectedOriginals;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder bytes=new ByteBuilder();
	private final ArrayList<Library> libraries=new ArrayList<Library>();
	private final ArrayList<Capture> captures=new ArrayList<Capture>();
	private static final int ALL_CONDITIONS=(1<<30)-1;
	private static final String[] SUMMARY_FIELDS={"full_library_class2_genes", "shard_class2_genes", "selected_saved_windows", "primary_model_alignments", "physical_edge_model_alignments"};
	private static final String BEST_HEADER="library\tcapture\toriginal_gene\tincomplete_coverage\tmissing44left20right\tstage\tcondition\tidentity\tmodel\tspan_bases\tforward_start0\tforward_stop0\tevent\tidentity_below064\tspan_below109\tfails_reciprocal_geometry\tgenuine_clip\tconstruction_model\tedge_identity_floor_pass\twinning_attempt_stage";
	private static final float[] CUTS={.60f,.62f,.64f};
	private static final String[] STAGES={"PRIMARY_RAW","PHYSICAL_EDGE_RAW","EITHER_RAW_NOT_FINAL_RECALL"};
	private static final String[] CONDITIONS={"BEST_POSITIVE_OVERLAP","OWN_MODEL_ANY_OVERLAP_STATE","OWN_MODEL_POSITIVE_OVERLAP",
		"ID060_SPAN109","ID060_NO_SPAN","ID062_SPAN109","ID062_NO_SPAN","ID064_SPAN109","ID064_NO_SPAN","OWN_ID064_SPAN109"};
	private static final String[] COHORTS={"CLASS2_ALL","CLASS2_INCOMPLETE","CLASS2_MISSING44LEFT20RIGHT"};
	private static final String[] QUANTITIES={"best_overlapping_identity","best_overlapping_span_bases","own_model_identity_anywhere","own_model_identity_overlapping"};
}
