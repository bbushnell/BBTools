package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.Outputs;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** Selects source panels and samples annotation-masked genomic windows.
 * Defaults balance genomes and draw eligible starts with replacement. Optional
 * lineage/GC balancing and nonoverlapping cores use a fixed length schedule;
 * exhaustion fails rather than changing that population. Coordinates are one-based
 * inclusive. One owner reads genomes and owns all output lifetimes. Source masking
 * never establishes a negative label: independent screening remains required.
 * @author Raiden
 */
public final class CobalaminShredSampler {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminShredSampler.class, false);
		try{parseAndRun(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses native flag=value options and validates the selected mode before output. */
	private static void parseAndRun(final String[] args){
		final CobalaminShredSampler x=new CobalaminShredSampler();
		String genomes=null, panel=null, panelIds=null, lengths=null, loci=null, out=null, groups=null;
		final Parser parser=new Parser();
		boolean threadsSet=false, overwriteSet=false, binsSet=false, panelCountSet=false, familiesSet=false;
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			require(b!=null && !b.isEmpty(), "Use flag=value: "+arg);
			if(a.equals("genomes")){
				genomes=b;
			}else if(a.equals("panel")){
				panel=b;
			}else if(a.equals("panelids")){
				panelIds=b;
			}else if(a.equals("lengths")){
				lengths=b;
			}else if(a.equals("loci")){
				loci=b;
			}else if(a.equals("outprefix")){
				out=b;
			}else if(a.equals("groups")){
				groups=b;
			}else if(a.equals("arm")){
				x.arm=b;
			}else if(a.equals("mode")){
				x.mode=b;
			}else if(a.equals("stratify")){
				x.stratify=b;
			}else if(a.equals("nonoverlap")){
				x.nonoverlap=Parse.parseBoolean(b);
			}else if(a.equals("gcbins")){
				x.gcBins=Parse.parseIntKMG(b);
				binsSet=true;
			}else if(a.equals("panelcount")){
				x.panelCount=Parse.parseIntKMG(b);
				panelCountSet=true;
			}else if(a.equals("count")){
				x.count=Parse.parseIntKMG(b);
			}else if(a.equals("flank")){
				x.flank=Parse.parseIntKMG(b);
			}else if(a.equals("genomecount")){
				x.genomeCount=Parse.parseIntKMG(b);
			}else if(a.equals("expectedgenomes")){
				x.expectedGenomes=Parse.parseIntKMG(b);
			}else if(a.equals("expectedpanel")){
				x.expectedPanel=Parse.parseIntKMG(b);
			}else if(a.equals("expectedlengths")){
				x.expectedLengths=Parse.parseIntKMG(b);
			}else if(a.equals("expectedloci")){
				x.expectedLoci=Parse.parseIntKMG(b);
			}else if(a.equals("seed")){
				x.seed=Long.parseLong(b);
			}else if(a.equals("families")){
				x.families=parseFamilies(b);
				familiesSet=true;
			}else if(a.equals("in") || a.equals("out") || a.equals("t") || a.equals("threads") || a.equals("ow") || a.equals("overwrite")){
				require(parser.parse(arg, a, b), "Unsupported standard option: "+arg);
				threadsSet|=a.equals("t") || a.equals("threads");
				overwriteSet|=a.equals("ow") || a.equals("overwrite");
			}else{throw new IllegalArgumentException("Unknown option: "+a);}
		}
		if(parser.in1!=null){
			require(lengths==null, "Specify in= or lengths= once");
			lengths=parser.in1;
		}
		if(parser.out1!=null){
			require(out==null, "Specify out= or outprefix= once");
			out=parser.out1;
		}
		require(!threadsSet || Shared.threads()==1, "Sampling has one owner; use t=1");
		require(!overwriteSet || !parser.overwrite, "Sampling preserves existing outputs; use a fresh prefix");
		require(x.mode.equals("sample") || x.mode.equals("panel"), "mode must be sample or panel");
		require(x.stratify.equals("none") || x.stratify.equals("phylumgc"), "stratify must be none or phylumgc");
		require(x.gcBins>0 && (!binsSet || x.stratify.equals("phylumgc")), "gcbins requires phylumgc and a positive bin count");
		Shared.setThreads(1);
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		if(x.mode.equals("panel")){
			require(genomes!=null && groups!=null && out!=null && x.expectedGenomes>0 && x.panelCount>0 && x.seed!=Long.MIN_VALUE,
				"Panel selection requires genomes, groups, outprefix, expectedgenomes, panelcount and seed");
			require(panel==null && panelIds==null && lengths==null && loci==null && x.arm==null && x.count<0 && x.flank<0
				&& x.genomeCount<0 && x.expectedPanel<0 && x.expectedLengths<0 && x.expectedLoci<0 && !x.nonoverlap && !familiesSet && x.stratify.equals("none"),
				"mode=panel selects a source roster; sampling options do not apply");
			x.selectPanel(genomes, groups, out);
			return;
		}
		require(groups==null && !panelCountSet, "groups/panelcount belong to mode=panel");
		require(genomes!=null && (panel!=null || panelIds!=null) && !(panel!=null && panelIds!=null) && lengths!=null && loci!=null && out!=null && x.arm!=null,
			"Bind genomes, panel= provenance ledger or panelids= source list, lengths, loci, outprefix and arm");
		require((x.arm.equals("panel") || x.arm.equals("nonpanel")) && x.count>0 && x.flank>=0 && x.genomeCount>0 && x.count>=x.genomeCount
			&& x.expectedGenomes>=x.genomeCount && x.expectedPanel>0 && x.expectedLengths>0 && x.expectedLoci>0 && x.seed!=Long.MIN_VALUE,
			"Declare positive population/sample counts, nonnegative flank, seed and panel/nonpanel arm");
		require(x.stratify.equals("none") || x.arm.equals("panel"), "phylumgc consumes an exact panel; use arm=panel");
		x.roster(genomes, panel, panelIds);
		x.lengths(lengths);
		x.masks(loci);
		if(x.stratify.equals("phylumgc")){x.prepareStrata();}
		x.write(out);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Selects whole, unambiguous source groups using only native metadata readers. */
	private void selectPanel(final String rosterPath, final String groupPath, final String prefix){
		final ObjectMap<String, PanelSource> sources=new ObjectMap<String, PanelSource>(String.class, PanelSource.class);
		final ObjectMap<String, PanelGroup> byGroup=new ObjectMap<String, PanelGroup>(String.class, PanelGroup.class);
		final ArrayList<PanelGroup> groups=new ArrayList<PanelGroup>();
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(groupPath)){
			RrnaResourceIO.header(input.reader, "genome\tassembly\ttaxid\tfasta_sha80\tsource_group\tpilot_exposed_group", groupPath);
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				require(p.terms()==6, "Source groups require the native six-field census schema");
				final String name=p.parseString(0), key=p.parseString(4);
				require(!name.isEmpty() && !key.isEmpty() && !sources.contains(name), "Source-group names must be nonempty and unique: "+name);
				PanelGroup group=byGroup.get(key);
				if(group==null){
					group=new PanelGroup(key);
					byGroup.put(key, group);
					groups.add(group);
				}
				group.members++;
				sources.put(name, new PanelSource(p.parseString(1), p.parseString(2), group));
			}
		}
		final ObjectSet<String> seen=new ObjectSet<String>(String.class), promised=new ObjectSet<String>(String.class);
		final ArrayList<String> lineages=new ArrayList<String>();
		try(Input input=new Input(rosterPath)){
			RrnaResourceIO.header(input.reader, NcrnaCorpusSample.GENOME_HEADER, rosterPath);
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				require(p.terms()==11, "Panel candidates require the complete original genome roster");
				final String name=p.parseString(1), lineage=p.parseString(4)+"/"+p.parseString(5);
				final PanelSource source=sources.get(name);
				require(source!=null && seen.add(name) && p.termEquals(source.assembly, 2) && p.termEquals(source.taxid, 3),
					"Panel source must join its census identity exactly once: "+name);
				if(promised.add(lineage)){lineages.add(lineage);}
				final PanelGroup group=source.group;
				group.present++;
				if(group.lineage==null){group.lineage=lineage;}
				else if(!group.lineage.equals(lineage)){group.ambiguous=true;}
				if(group.representative==null || name.compareTo(group.representative)<0){
					group.representative=name;
					group.row=row;
				}
			}
		}
		require(seen.size()==expectedGenomes, "Panel candidate genome count changed: "+seen.size()+" expected="+expectedGenomes);
		final ObjectMap<String, PanelStratum> byLineage=new ObjectMap<String, PanelStratum>(String.class, PanelStratum.class);
		final ArrayList<PanelStratum> strata=new ArrayList<PanelStratum>();
		int eligible=0;
		for(PanelGroup group:groups){
			if(group.present==0){continue;}
			require(group.present==group.members, "Source groups must enter or leave the candidate roster whole: "+group.key);
			if(group.ambiguous){continue;}
			PanelStratum stratum=byLineage.get(group.lineage);
			if(stratum==null){
				stratum=new PanelStratum(group.lineage);
				byLineage.put(group.lineage, stratum);
				strata.add(stratum);
			}
			stratum.groups.add(group);
			eligible++;
		}
		for(String lineage:lineages){
			require(byLineage.contains(lineage), "No unambiguous source group can represent promised phylum: "+lineage);
		}
		Collections.sort(strata, (a, b)->a.key.compareTo(b.key));
		final int target=Math.min(panelCount, eligible), remaining=target-strata.size(), capacity=eligible-strata.size();
		require(target>0 && remaining>=0, "Panel cap cannot cover every promised phylum: cap="+panelCount+" phyla="+strata.size());
		int assigned=0;
		for(PanelStratum stratum:strata){
			final long weight=(long)remaining*(stratum.groups.size()-1);
			stratum.quota=1+(capacity==0 ? 0 : (int)(weight/capacity));
			stratum.remainder=capacity==0 ? 0 : weight%capacity;
			assigned+=stratum.quota;
		}
		final ArrayList<PanelStratum> remainderOrder=new ArrayList<PanelStratum>(strata);
		Collections.sort(remainderOrder, (a, b)->a.remainder==b.remainder ? a.key.compareTo(b.key) : Long.compare(b.remainder, a.remainder));
		for(PanelStratum stratum:remainderOrder){
			if(assigned<target && stratum.quota<stratum.groups.size()){
				stratum.quota++;
				assigned++;
			}
		}
		require(assigned==target, "Capacity-weighted panel quotas must conserve the requested slots");
		final Random random=new Random(seed);
		for(PanelStratum stratum:strata){
			Collections.sort(stratum.groups, (a, b)->a.key.compareTo(b.key));
			Collections.shuffle(stratum.groups, random);
			for(int i=0; i<stratum.quota; i++){stratum.groups.get(i).selected=true;}
		}
		Collections.sort(groups, (a, b)->a.key.compareTo(b.key));
		writePanel(prefix, groups, strata, target, eligible);
	}

	/** Writes the exact panel and every candidate group's disposition before genome I/O. */
	private void writePanel(final String prefix, final ArrayList<PanelGroup> groups, final ArrayList<PanelStratum> strata,
			final int target, final int eligible){
		assert(target>0 && target<=eligible) : "Panel output must retain the selected subset of eligible groups";
		RrnaResourceIO.fresh(prefix, ".ids.txt", ".genomes.tsv", ".groups.tsv", ".summary.tsv");
		try(Outputs owners=new Outputs(4)){
			final ByteStreamWriter ids=textOutput(owners, prefix+".ids.txt"), genomes=textOutput(owners, prefix+".genomes.tsv");
			final ByteStreamWriter rows=textOutput(owners, prefix+".groups.tsv"), summary=textOutput(owners, prefix+".summary.tsv");
			final ByteBuilder b=new ByteBuilder();
			genomes.println(NcrnaCorpusSample.GENOME_HEADER);
			rows.println("source_group\tstratum\tsources\trepresentative\tselected\treason");
			int written=0;
			for(PanelGroup group:groups){
				if(group.present==0){continue;}
				if(group.selected){
					ids.println(group.representative);
					genomes.println(group.row);
					written++;
				}
				final String reason=group.ambiguous ? "AMBIGUOUS_TAXONOMY" : group.selected ? "SELECTED" : "NOT_SELECTED";
				rows.print(b.clear().append(group.key).tab().append(group.ambiguous ? "MULTIPLE" : group.lineage).tab().append(group.present)
					.tab().append(group.representative).tab().append(group.selected).tab().append(reason).nl());
			}
			require(written==target, "Panel IDs and original roster rows must conserve selected groups");
			summary.println("candidate_genomes\teligible_groups\tphyla\tpanel_groups\tseed\tlabel");
			summary.print(b.clear().append(expectedGenomes).tab().append(eligible).tab().append(strata.size()).tab().append(written)
				.tab().append(seed).append("\tPROVISIONAL_SOURCE_PANEL\n"));
		}
	}

	/** Measures exact source-genome GC and allocates fixed strata/source quotas. */
	private void prepareStrata(){
		assert(selected.size()==expectedPanel) : "Stratified sampling consumes the exact supplied panel";
		Collections.sort(selected, (a, b)->a.name.compareTo(b.name));
		final ObjectMap<String, BackgroundStratum> byStratum=new ObjectMap<String, BackgroundStratum>(String.class, BackgroundStratum.class);
		for(Genome genome:selected){
			final long[] counts=new long[2];
			final ObjectSet<String> names=new ObjectSet<String>(String.class);
			RrnaResourceIO.stream(genome.fasta, read->{
				final String name=RrnaResourceIO.first(read.id);
				require(read.bases!=null && names.add(name), "GC census needs unique contigs: "+genome.name+"/"+name);
				final Mask mask=genome.masks.get(name);
				require(mask==null || mask.length==read.length(), "Source-mask length disagrees with GC-pass FASTA: "+genome.name+"/"+name);
				for(byte base:read.bases){
					final int code=RrnaResourceIO.base(base);
					if(code>=0){counts[0]++;}
					if(code==1 || code==2){counts[1]++;}
				}
			});
			for(String name:genome.maskNames){require(names.contains(name), "Source mask contig missing during GC census: "+genome.name+"/"+name); }
			require(counts[0]>0, "Cannot assign a GC stratum without defined bases: "+genome.name);
			genome.defined=counts[0];
			genome.gc=counts[1];
			genome.gcBin=gcBin(genome.gc, genome.defined, gcBins);
			final String key=genome.lineage+"/gc"+genome.gcBin;
			BackgroundStratum stratum=byStratum.get(key);
			if(stratum==null){
				stratum=new BackgroundStratum(key);
				byStratum.put(key, stratum);
				backgroundStrata.add(stratum);
			}
			stratum.genomes.add(genome);
			genome.stratum=stratum;
		}
		Collections.sort(backgroundStrata, (a, b)->a.key.compareTo(b.key));
		require(count>=backgroundStrata.size(), "Window count cannot cover every occupied panel stratum");
		final Random allocation=new Random(seed^0x391de52bL);
		final ArrayList<BackgroundStratum> extra=new ArrayList<BackgroundStratum>(backgroundStrata);
		Collections.shuffle(extra, allocation);
		for(BackgroundStratum stratum:backgroundStrata){stratum.quota=count/backgroundStrata.size(); }
		for(int i=0; i<count%backgroundStrata.size(); i++){extra.get(i).quota++;}
		for(BackgroundStratum stratum:backgroundStrata){
			final ArrayList<Genome> sourceOrder=new ArrayList<Genome>(stratum.genomes);
			Collections.shuffle(sourceOrder, allocation);
			for(Genome genome:sourceOrder){genome.quota=stratum.quota/sourceOrder.size(); }
			for(int i=0; i<stratum.quota%sourceOrder.size(); i++){sourceOrder.get(i).quota++;}
		}
		int offset=0;
		for(Genome genome:selected){
			genome.scheduleStart=offset;
			offset+=genome.quota;
		}
		require(offset==count, "Stratum/source allocation must consume every fixed length-schedule entry");
	}

	/** Exact half-open composition bin; the all-GC endpoint belongs to the final bin. */
	public static int gcBin(final long gc, final long defined, final int bins){
		require(defined>0 && gc>=0 && gc<=defined && bins>0, "GC bins require defined bases, bounded GC count and positive bins");
		return gc==defined ? bins-1 : (int)(Math.multiplyExact(gc, (long)bins)/defined);
	}

	/** Owns a fresh text writer using the legacy backend flags. */
	private static ByteStreamWriter textOutput(final Outputs owners, final String path){
		return owners.open(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));
	}

	/** Resolves exact panel membership; the legacy nonpanel arm excludes assembly/taxid aliases. */
	private void roster(final String path, final String panelPath, final String panelIdsPath){
		final ObjectSet<String> panelFiles=new ObjectSet<String>(String.class), panelAssemblies=new ObjectSet<String>(String.class), panelTaxids=new ObjectSet<String>(String.class);
		final ObjectSet<String> panelNames=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(FileFormat.testInput(panelPath==null ? panelIdsPath : panelPath, FileFormat.TEXT, null, true, true))){
			final ByteFile in=input.reader;
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);
				if(panelPath==null){
					require(p.terms()==1 && row.length>0 && panelNames.add(p.parseString(0)), "Panel IDs must be unique source names");
				}else{
					require(p.terms()==2 && p.parseString(0).matches("[0-9a-f]{20}") && p.parseString(1).startsWith("/") && panelFiles.add(p.parseString(1)),
						"Panel paths must be unique absolute FASTA identities with sha80 provenance");
				}
			}
		}
		require((panelPath==null ? panelNames.size() : panelFiles.size())==expectedPanel, "Panel population changed");
		final ArrayList<Genome> all=new ArrayList<Genome>();
		final ObjectSet<String> seen=new ObjectSet<String>(String.class), paths=new ObjectSet<String>(String.class);
		int matched=0;
		try(Input input=new Input(FileFormat.testInput(path, FileFormat.TEXT, null, true, true))){
			final ByteFile in=input.reader;
			RrnaResourceIO.header(in, NcrnaCorpusSample.GENOME_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);
				require(p.terms()==11, "Expected complete source genome roster schema");
				final Genome g=new Genome(p.parseString(1), p.parseString(2), p.parseString(3), p.parseString(7));
				g.lineage=p.parseString(4)+"/"+p.parseString(5);
				require(g.fasta.startsWith("/") && seen.add(g.name) && paths.add(g.fasta), "Unique genome names and absolute FASTA paths are required");
				g.panel=panelPath==null ? panelNames.contains(g.name) : panelFiles.contains(g.fasta);
				all.add(g);
				if(g.panel){
					matched++;
					panelAssemblies.add(stem(g.assembly));
					panelTaxids.add(g.taxid);
				}
			}
		}
		require(all.size()==expectedGenomes && matched==expectedPanel, "Full roster and panel membership must match declared populations");
		for(Genome g:all){
			if(arm.equals("panel")?g.panel:!g.panel && !panelAssemblies.contains(stem(g.assembly)) && !panelTaxids.contains(g.taxid)){selected.add(g);}
		}
		require(selected.size()>=genomeCount && (!arm.equals("panel") || genomeCount==expectedPanel), "Panel arm uses every panel genome; nonpanel arm needs enough source-disjoint genomes");
		eligibleGenomes=selected.size();
		Collections.sort(selected, (a, b)->a.name.compareTo(b.name));
		Collections.shuffle(selected, new Random(seed^0x6732a74bL));
		while(selected.size()>genomeCount){selected.remove(selected.size()-1);}
		for(Genome g:selected){byGenome.put(g.name, g); }
	}
	/** Draws once from sorted empirical lengths; downstream failure never redraws. */
	private void lengths(final String path){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(path, r->{
			require(r.bases!=null && r.length()>0 && ids.add(RrnaResourceIO.first(r.id)), "Empirical lengths require unique intact positive records");
			lengths.add(r.length());
		});
		require(lengths.size==expectedLengths, "Empirical positive population changed: "+lengths.size+" expected="+expectedLengths);
		// Sorting makes schedules independent of FASTA row order without changing frequencies.
		Arrays.sort(lengths.array, 0, lengths.size);
		final Random random=new Random(seed);
		schedule=new int[count];
		for(int i=0; i<count; i++){schedule[i]=lengths.get(random.nextInt(lengths.size));}
	}
	/** Validates the complete locus population and retains masks only for selected sources. */
	private void masks(final String files){
		final String[] paths=files.split(",");
		require(paths.length==families.length, "Supply one source locus ledger per declared family, in the same order");
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		int observed=0;
		for(int f=0; f<families.length; f++){
			try(Input input=new Input(FileFormat.testInput(paths[f], FileFormat.TEXT, null, true, true))){
				final ByteFile in=input.reader;
				RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER, paths[f]);
				for(byte[] row; (row=in.nextLine())!=null; ){
					p.set(row);
					require(p.terms()==21, "All source masks need the original annotation schema");
					final String id=p.parseString(0);
					require(id.startsWith(families[f]+"_") && ids.add(id), "Source mask IDs must retain unique declared family provenance");
					observed++;
					final Genome g=byGenome.get(p.parseString(1));
					if(g==null){continue;}
					require(p.termEquals(g.assembly, 2) && p.termEquals(g.taxid, 3), "Mask source identity must match selected genome: "+id);
					final String contig=p.parseString(6);
					Mask mask=g.masks.get(contig);
					if(mask==null){
						mask=new Mask(p.parseInt(15));
						g.masks.put(contig, mask);
						g.maskNames.add(contig);
					}
					require(mask.length==p.parseInt(15) && mask.length>0, "Source masks must agree on contig length: "+id);
					// Mask the full annotation span, including any internal interval gaps.
					int lo=Integer.MAX_VALUE, hi=-1;
					for(String interval:p.parseString(8).split(",")){
						final int dash=interval.indexOf('-');
						require(dash>0, "Source interval must contain inclusive bounds: "+id);
						final int a=Integer.parseInt(interval.substring(0, dash)), z=Integer.parseInt(interval.substring(dash+1));
						require(a>0 && z>=a && z<=mask.length, "Source mask must lie within its contig: "+id);
						lo=Math.min(lo, a);
						hi=Math.max(hi, z);
					}
					mask.spans.add(new int[]{lo, hi});
					g.sourceLoci++;
				}
			}
		}
		require(observed==expectedLoci, "The complete declared-family annotation mask population changed");
	}
	/** An explicit family list changes masks only; legacy three-family defaults remain exact. */
	private static String[] parseFamilies(final String value){
		assert(value!=null):"The families= option supplies an explicit ordered list matching loci=";
		final String[] names=value.split(",", -1);
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		for(String name:names){require(name.matches("RF[0-9]{5}") && seen.add(name), "Family IDs must be unique Rfam accessions in source-ledger order: "+name); }
		return names;
	}
	/** Writes one fixed schedule; every opened output is owned even on partial failure. */
	private void write(String prefix){
		assert(schedule!=null && schedule.length==count && selected.size()==genomeCount) :
			"Window output requires the complete fixed schedule and resolved source allocation";
		final boolean stratified=stratify.equals("phylumgc"), extended=stratified || nonoverlap;
		RrnaResourceIO.fresh(prefix, ".cores.fa.gz", ".contexts.fa.gz", ".shreds.tsv", ".genomes.tsv", ".lengths.tsv", ".summary.tsv");
		if(extended){
			RrnaResourceIO.fresh(prefix, ".schedule.tsv", ".window_metadata.tsv");
			if(stratified){RrnaResourceIO.fresh(prefix, ".strata.tsv");}
			writeSchedule(prefix, stratified);
		}
		final ByteBuilder b=new ByteBuilder();
		final Random positions=new Random(seed^0x1b93f045L);
		final Selection chosen=new Selection();
		long coreBases=0, contextBases=0;
		int written=0;
		try(Outputs owners=new Outputs(6);
			SequenceOutput core=new SequenceOutput(prefix+".cores.fa.gz", 0);
			SequenceOutput context=new SequenceOutput(prefix+".contexts.fa.gz", flank)){
			final ByteStreamWriter rows=textOutput(owners, prefix+".shreds.tsv"), genomes=textOutput(owners, prefix+".genomes.tsv");
			final ByteStreamWriter hist=textOutput(owners, prefix+".lengths.tsv"), summary=textOutput(owners, prefix+".summary.tsv");
			final ByteStreamWriter metadata=extended ? textOutput(owners, prefix+".window_metadata.tsv") : null;
			final ByteStreamWriter strata=stratified ? textOutput(owners, prefix+".strata.tsv") : null;
			if(metadata!=null){metadata.println("shred\tgenome\tstratum\tgenome_gc_bin\tgenome_defined_bases\tgenome_gc_bases\twindow_gc_bases\twindow_bases");}
			rows.println(HEADER);
			genomes.println("arm\tgenome\tassembly\ttaxid\tfasta\tsource_mask_loci\tshreds\tbases");
			for(int gi=0; gi<selected.size(); gi++){
				final Genome g=selected.get(gi);
				final ArrayList<Contig> contigs=new ArrayList<Contig>();
				final ObjectSet<String> names=new ObjectSet<String>(String.class);
				if(stratified && g.quota==0){
					genomes.print(b.clear().append(arm).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(g.fasta)
						.tab().append(g.sourceLoci).append("\t0\t0\n"));
					continue;
				}
				RrnaResourceIO.stream(g.fasta, r->{
					final String name=RrnaResourceIO.first(r.id);
					require(r.bases!=null && names.add(name), "Source FASTA needs unique named contigs");
					final Mask source=g.masks.get(name);
					require(source==null || source.length==r.length(), "Mask coordinates must match actual FASTA contig length: "+name);
					final Euk18sShredSet.Contig mask=new Euk18sShredSet.Contig(name, r.length(), 0, 1, 1);
					if(source!=null){mask.mask.addAll(source.spans);}
					for(int i=0; i<r.length(); ){
						if(RrnaResourceIO.base(r.bases[i])>=0){
							i++;
							continue;
						}
						final int start=i++;
						while(i<r.length() && RrnaResourceIO.base(r.bases[i])<0){i++;}
						mask.mask.add(new int[]{start+1, i});
					}
					mask.prepareMask(0);
					contigs.add(new Contig(name, r.bases, mask));
				});
				if(extended){Collections.sort(contigs, (a, z)->a.name.compareTo(z.name));}
				for(String name:g.maskNames){require(names.contains(name), "Every selected source annotation contig must occur in its FASTA: "+name); }
				int n=0;
				long bases=0;
				final int first=stratified ? g.scheduleStart : gi, limit=stratified ? g.scheduleStart+g.quota : count;
				final int step=stratified ? 1 : genomeCount;
				for(int i=first; i<limit; i+=step){
					final int length=schedule[i];
					final Contig c;
					final int lo;
					if(extended){
						final long total=walkStarts(contigs, length, -1, null);
						require(total>0, "Fixed schedule exhausted: window="+i+" source="+g.name+" stratum="+stratumName(g)+" length="+length);
						chosen.contig=null;
						walkStarts(contigs, length, Euk18sShredSet.bounded(positions, total), chosen);
						require(chosen.contig!=null, "An unchanged eligible-start universe must contain its drawn ticket");
						c=chosen.contig;
						lo=chosen.start;
					}else{
						// Preserve the legacy path's exact traversal and random-number consumption.
						final ArrayList<Start> starts=startSpace(contigs, length);
						require(!starts.isEmpty(), "Assigned length has no unmasked start in genome; do not redraw its length: "+g.name+" length="+length);
						final long ticket=Euk18sShredSet.bounded(positions, starts.get(starts.size()-1).end);
						int a=0, z=starts.size();
						while(a<z){
							final int mid=(a+z)>>>1;
							if(starts.get(mid).end<=ticket){a=mid+1;}else{z=mid;}
						}
						final Start s=starts.get(a);
						c=s.contig;
						lo=s.start+(int)(ticket-(a==0 ? 0 : starts.get(a-1).end));
					}
					final int hi=lo+length-1;
					require(!c.mask.overlaps(lo, hi), "Complement-start sampler must exclude all source annotations and ambiguous core bases");
					if(nonoverlap){addInterval(c.mask, lo, hi);}
					final int left=(int)Math.max(1L, (long)lo-flank), right=(int)Math.min(c.bases.length, (long)hi+flank);
					final String id=arm+"_"+i;
					core.add(new Read(Arrays.copyOfRange(c.bases, lo-1, hi), null, id, i));
					context.add(new Read(Arrays.copyOfRange(c.bases, left-1, right), null, id, i));
					rows.print(b.clear().append(id).tab().append(arm).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(c.name)
						.tab().append(lo).tab().append(hi).tab().append(left).tab().append(right).tab().append(lo-left+1).tab().append(hi-left+1).tab().append(length).nl());
					if(metadata!=null){
						int gc=0;
						for(int pos=lo-1; pos<hi; pos++){
							final int code=RrnaResourceIO.base(c.bases[pos]);
							if(code==1 || code==2){gc++;}
						}
						b.clear().append(id).tab().append(g.name).tab().append(stratumName(g)).tab();
						if(stratified){b.append(g.gcBin).tab().append(g.defined).tab().append(g.gc);}
						else{b.append("NA\tNA\tNA");}
						metadata.print(b.tab().append(gc).tab().append(length).nl());
					}
					n++;
					written++;
					bases+=length;
					coreBases+=length;
					contextBases+=right-left+1;
				}
				if(stratified){
					require(n==g.quota, "Every source must realize its fixed window quota: "+g.name);
					g.stratum.written+=n;
					if(n>0){g.stratum.contributors++;}
				}
				genomes.print(b.clear().append(arm).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(g.fasta).tab().append(g.sourceLoci).tab().append(n).tab().append(bases).nl());
			}
			require(written==count, "Every fixed length schedule entry must become one shred");
			final IntList empirical=new IntList(), sampled=new IntList();
			for(int i=0; i<lengths.size; i++){empirical.increment(lengths.get(i), 1);}
			for(int n:schedule){sampled.increment(n, 1); }
			hist.println("length\tpositive_corpus\tbackground_sample");
			for(int i=0; i<Math.max(empirical.size, sampled.size); i++){
				if(empirical.get(i)>0 || sampled.get(i)>0){hist.print(b.clear().append(i).tab().append(empirical.get(i)).tab().append(sampled.get(i)).nl());}
			}
			summary.println("arm\teligible_genomes\tsampled_genomes\tshreds\tcore_bases\tcontext_bases\tseed\tflank\tpolicy\tlabel");
			b.clear().append(arm).tab().append(eligibleGenomes).tab().append(genomeCount).tab().append(written).tab().append(coreBases).tab().append(contextBases).tab().append(seed).tab().append(flank);
			if(!extended){b.append("\tBALANCED_GENOMES_EMPIRICAL_LENGTH_UNIFORM_ELIGIBLE_START_WITH_REPLACEMENT");}
			else{b.tab().append(stratified ? "BALANCED_PANEL_PHYLA_GC" : "BALANCED_GENOMES").append(nonoverlap ? "_NONOVERLAPPING_FIXED_LENGTHS" : "_FIXED_LENGTHS_WITH_REPLACEMENT");}
			summary.print(b.append("\tSOURCE_MASKED_CM_SCREEN_PENDING\n"));
			if(strata!=null){
				strata.println("stratum\tselected_sources\tcontributing_sources\tscheduled_windows\twritten_windows");
				for(BackgroundStratum stratum:backgroundStrata){
					require(stratum.written==stratum.quota, "Occupied-stratum quotas must reconcile: "+stratum.key);
					strata.print(b.clear().append(stratum.key).tab().append(stratum.genomes.size()).tab().append(stratum.contributors)
						.tab().append(stratum.quota).tab().append(stratum.written).nl());
				}
			}
		}
		System.err.println("COBALAMIN_SHRED_SAMPLE_PASS arm="+arm+" shreds="+written+" bases="+coreBases+"; all candidates still require both-strand clan CM screening");
	}

	/** Saves the complete target schedule before any window can fail to materialize. */
	private void writeSchedule(final String prefix, final boolean stratified){
		assert(schedule!=null && schedule.length==count) : "The exhaustion record must include every requested window";
		try(Output owner=new Output(prefix+".schedule.tsv")){
			final ByteBuilder b=new ByteBuilder();
			owner.writer.println("window\tgenome\tstratum\tlength");
			int source=0;
			for(int i=0; i<count; i++){
				if(stratified){
					while(source<selected.size() && i>=selected.get(source).scheduleStart+selected.get(source).quota){source++;}
					require(source<selected.size(), "Every scheduled length has a source owner");
				}
				final Genome genome=selected.get(stratified ? source : i%genomeCount);
				owner.writer.print(b.clear().append(arm).append('_').append(i).tab().append(genome.name).tab().append(stratumName(genome)).tab().append(schedule[i]).nl());
			}
		}
	}

	private static String stratumName(final Genome genome){return genome.stratum==null ? "ALL" : genome.stratum.key;}

	/** Counts eligible starts, or stops after locating a ticket without allocating gap objects. */
	private static long walkStarts(final ArrayList<Contig> contigs, final int length, final long ticket, final Selection selected){
		assert(length>0 && (selected==null || ticket>=0)) : "Start-space walks use a positive scheduled length and a drawn nonnegative ticket";
		long total=0;
		for(Contig contig:contigs){
			final long last=(long)contig.bases.length-length+1;
			if(last<1){continue;}
			long next=1;
			for(int i=0; i<=contig.mask.starts.size; i++){
				final boolean tail=i==contig.mask.starts.size;
				final long from=tail ? last+1 : Math.max(1L, (long)contig.mask.starts.get(i)-length+1);
				final long end=Math.min(last, from-1);
				if(end>=next){
					final long after=Math.addExact(total, end-next+1);
					if(selected!=null && ticket<after){
						selected.contig=contig;
						selected.start=(int)(next+ticket-total);
						return after;
					}
					total=after;
				}
				if(tail){break;}
				next=Math.max(next, (long)contig.mask.ends.get(i)+1);
				if(next>last){break;}
			}
		}
		return total;
	}

	/** Adds an inclusive exclusion to the prepared primitive union; never reinitializes it. */
	public static void addInterval(final Euk18sShredSet.Contig mask, int lo, int hi){
		require(lo>0 && hi>=lo && hi<=mask.length, "Incremental exclusion must lie inside its contig: "+lo+"-"+hi);
		final IntList starts=mask.starts, ends=mask.ends;
		assert(starts.size==ends.size) : "Start-space enumeration requires paired sorted interval bounds";
		final int size=starts.size;
		int first=0;
		while(first<size && (long)ends.get(first)+1<lo){first++;}
		int after=first;
		while(after<size && starts.get(after)<=(long)hi+1){
			lo=Math.min(lo, starts.get(after));
			hi=Math.max(hi, ends.get(after));
			after++;
		}
		final int newSize=size-(after-first)+1;
		starts.ensureCapacity(Math.max(0, newSize-size));
		ends.ensureCapacity(Math.max(0, newSize-size));
		System.arraycopy(starts.array, after, starts.array, first+1, size-after);
		System.arraycopy(ends.array, after, ends.array, first+1, size-after);
		starts.array[first]=lo;
		ends.array[first]=hi;
		starts.setSize(newSize);
		ends.setSize(newSize);
		assert(mask.overlaps(lo, hi)) : "The next scheduled length must see each newly excluded core";
	}
	/** Retains the historical cumulative-gap traversal for byte-identical default sampling. */
	private static ArrayList<Start> startSpace(final ArrayList<Contig> contigs, final int length){
		require(length>0, "Complement start-space requires positive fixed core length");
		final ArrayList<Start> result=new ArrayList<Start>();
		long total=0;
		for(Contig c:contigs){
			final long last=(long)c.bases.length-length+1;
			if(last<1){continue;}
			long next=1;
			for(int i=0; i<c.mask.starts.size; i++){
				final long from=Math.max(1L, (long)c.mask.starts.get(i)-length+1), to=Math.min(last, c.mask.ends.get(i));
				if(from>next){
					final long end=Math.min(last, from-1);
					if(end>=next){
						total=Math.addExact(total, end-next+1);
						result.add(new Start(c, (int)next, total));
					}
				}
				next=Math.max(next, to+1);
				if(next>last){break;}
			}
			if(next<=last){
				total=Math.addExact(total, last-next+1);
				result.add(new Start(c, (int)next, total));
			}
		}
		return result;
	}

	/** Assembly-version aliases remain excluded in the unchanged legacy nonpanel arm. */
	private static String stem(String assembly){
		require(!assembly.isEmpty(), "Assembly provenance cannot be empty");
		final int dot=assembly.indexOf('.');
		return dot<0 ? assembly : assembly.substring(0, dot);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	/** Selected identity, annotation masks and optional fixed stratum allocation. */
	private static final class Genome {
		Genome(String n, String a, String t, String f){
			name=n;
			assembly=a;
			taxid=t;
			fasta=f;
		}
		final String name, assembly, taxid, fasta;
		String lineage;
		boolean panel;
		int sourceLoci, gcBin, quota, scheduleStart;
		long defined, gc;
		BackgroundStratum stratum;
		final ObjectMap<String, Mask> masks=new ObjectMap<String, Mask>(String.class, Mask.class);
		final ArrayList<String> maskNames=new ArrayList<String>();
	}

	/** Owns the existing one-thread FASTA writer; resource scopes preserve body failures. */
	private static final class SequenceOutput implements AutoCloseable {
		SequenceOutput(final String path, final int flank){
			assert(Shared.threads()==1) : "The existing FASTA output path is direct and unordered for this single-owner sampler";
			output=new Euk18sShredSet.Output(path, flank);
		}
		void add(final Read read){output.add(read);}
		@Override
		public void close(){
			if(output.close()){throw new IllegalStateException("Background FASTA write failed");}
		}
		final Euk18sShredSet.Output output;
	}

	/** Census identity and its indivisible source group, not biological truth. */
	private static final class PanelSource {
		PanelSource(final String a, final String t, final PanelGroup g){
			assembly=a;
			taxid=t;
			group=g;
		}
		final String assembly, taxid;
		final PanelGroup group;
	}

	/** Lexical representative row and membership/disposition counters for one group. */
	private static final class PanelGroup {
		PanelGroup(final String name){key=name;}
		final String key;
		String representative, lineage;
		byte[] row;
		int members, present;
		boolean ambiguous, selected;
	}

	/** Eligible groups in one domain/phylum, with largest-remainder allocation state. */
	private static final class PanelStratum {
		PanelStratum(final String name){key=name;}
		final String key;
		final ArrayList<PanelGroup> groups=new ArrayList<PanelGroup>();
		int quota;
		long remainder;
	}

	/** Sources sharing domain/phylum/GC bin; contributors may be fewer than sources. */
	private static final class BackgroundStratum {
		BackgroundStratum(final String name){key=name;}
		final String key;
		final ArrayList<Genome> genomes=new ArrayList<Genome>();
		int quota, written, contributors;
	}

	/** Reused result slot for allocation-free start-space lookup in the new modes. */
	private static final class Selection {
		Contig contig;
		int start;
	}
	/** Inclusive source annotations retained before loading the matching contig. */
	private static final class Mask {
		Mask(int n){length=n;}
		final int length;
		final ArrayList<int[]> spans=new ArrayList<int[]>();
	}
	/** Owns a loaded sequence and its prepared exclusion union for one source genome. */
	private static final class Contig {
		Contig(String n, byte[] b, Euk18sShredSet.Contig m){
			name=n;
			bases=b;
			mask=m;
		}
		final String name;
		final byte[] bases;
		final Euk18sShredSet.Contig mask;
	}
	/** Historical start span; end is the exclusive cumulative ticket bound. */
	private static final class Start {
		Start(Contig c, int s, long e){
			contig=c;
			start=s;
			end=e;
		}
		final Contig contig;
		final int start;
		final long end;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	private final ArrayList<Genome> selected=new ArrayList<Genome>();
	private final ArrayList<BackgroundStratum> backgroundStrata=new ArrayList<BackgroundStratum>();
	private final ObjectMap<String, Genome> byGenome=new ObjectMap<String, Genome>(String.class, Genome.class);
	private final IntList lengths=new IntList();
	/** Immutable requested lengths indexed by emitted window ID. */
	private int[] schedule;
	private String arm;
	private String mode="sample", stratify="none";
	private boolean nonoverlap=false;
	private int gcBins=10, panelCount=1000;
	private int count=-1, flank=-1, genomeCount=-1, expectedGenomes=-1, expectedPanel=-1, expectedLengths=-1, expectedLoci=-1, eligibleGenomes;
	private long seed=Long.MIN_VALUE;
	private static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	private String[] families=FAMILIES;
	public static final String HEADER="shred\tarm\tgenome\tassembly\ttaxid\tseqid\tcore_start1\tcore_stop1\tcontext_start1\tcontext_stop1\tcore_in_context_start1\tcore_in_context_stop1\tlength";
}
