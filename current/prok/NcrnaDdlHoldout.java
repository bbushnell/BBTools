package prok;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;

/** Current-roster assembly/alias and exact-literal-DNA exclusions for DDL votes.
 * File identities come from the formal roster and reviewed alias groups, never
 * TaxID. All panel genes are observed before any retained member is admitted.
 * The sample additionally omits genomes containing an excluded duplicate voter.
 * Metadata is loaded once; sequence arrays are borrowed until voting finishes.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlHoldout implements NcrnaDdlMembers.RosterFilter{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlHoldout(final String manifest, final String panel, final String sources, final String aliases){
		this(manifest, panel, sources, aliases, "RF00174");
	}

	/** Binds observations to one actual clan member, preserving the old main-member API. */
	public NcrnaDdlHoldout(final String manifest, final String panel, final String sources,
			final String aliases, final String member_){
		this(manifest, panel, sources, aliases, member_, null, 0);
	}

	/** Appends reviewed historical sources without inventing native source-ledger IDs. */
	public NcrnaDdlHoldout(final String manifest, final String panel, final String sources,
			final String aliases, final String member_, final String protectedPriors, final int expectedPriors){
		this(manifest, panel, sources, aliases, member_, protectedPriors, expectedPriors, false);
	}

	/** Release refitting keeps the checked source join but admits every eligible voter.
	 * No evaluation panel is read; aliases and literal sequence equality exclude none. */
	public NcrnaDdlHoldout(final String manifest, final String panel, final String sources,
			final String aliases, final String member_, final String protectedPriors,
			final int expectedPriors, final boolean refit_){
		require(!refit_ || panel==null, "All-data refitting omits panel=; it cannot silently ignore an evaluation panel");
		refit=refit_;
		require(protectedPriors==null ? expectedPriors==0 : expectedPriors>0,
			"Protected holdout sources require the companion and its exact positive count together");
		require("RF00174".equals(member_) || "RF01482".equals(member_) || "RF01689".equals(member_),
			"Holdout sources require an explicit supported cobalamin member identity");
		member=member_;
		final ObjectIntMap<String> groups=readAliases(aliases);
		readGenomes(manifest, groups);
		parent=new int[genomes.size()];
		for(int i=0; i<parent.length; i++){parent[i]=i;}
		final ObjectIntMap<String> tokens=new ObjectIntMap<String>(String.class);
		for(int i=0; i<genomes.size(); i++){
			final Genome g=genomes.get(i);
			connect(tokens, g.canonical, i);
			if(present(g.accession)){connect(tokens, "accession:"+g.accession, i);}
			if(g.alias>=0){connect(tokens, "alias:"+g.alias, i);}
		}
		held=new boolean[parent.length];
		duplicateGenome=new boolean[parent.length];
		sampled=new boolean[parent.length];
		if(!refit){readPanel(panel);}
		readSources(sources);
		if(protectedPriors!=null){readProtected(protectedPriors, expectedPriors);}
	}

	private static ObjectIntMap<String> readAliases(final String path){
		final ObjectIntMap<String> result=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GROUP_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==14 && p.parseInt(10)>=0, "Aliases need complete reviewed resolved-group rows");
				for(int column=2; column<=3; column++){
					final String key=p.parseString(column);
					final int group=p.parseInt(10);
					require(!result.contains(key) || result.get(key)==group, "One literal/canonical FASTA cannot change alias groups");
					result.put(key, group);
				}
			}
		}
		return result;
	}

	private void readGenomes(final String path, final ObjectIntMap<String> aliases){
		NcrnaLearningExclusions.readManifest(path, byAssembly, new ObjectIntMap<String>(String.class),
			new ArrayList<String>(), (id, fasta, source, row)->{
				final String canonical=canonical(fasta);
				int group=aliases.get(source);
				for(final String pathAlias:new String[]{fasta, canonical}){
					final int other=aliases.get(pathAlias);
					require(group<0 || other<0 || group==other,
						"Source, linked and canonical FASTA aliases have conflicting reviewed groups: "+id);
					if(group<0){group=other;}
				}
				if(group<0){unmappedAliases++;}
				genomes.add(new Genome(id, row.parseString(1), fasta, row.parseString(4),
					row.parseString(6), canonical, group));
			});
		require(!genomes.isEmpty(), "DDL preparation needs the explicit nonempty formal genome population");
	}

	private void connect(final ObjectIntMap<String> tokens, final String token, final int genome){
		final int previous=tokens.get(token);
		if(previous<0){tokens.put(token, genome);}
		else{
			final int a=root(previous), b=root(genome);
			if(a!=b){parent[b]=a;}
		}
	}

	private int root(int index){
		assert(index>=0 && index<parent.length) : "Assembly aliases refer only to bound formal-roster ordinals";
		while(parent[index]!=index){
			parent[index]=parent[parent[index]];
			index=parent[index];
		}
		return index;
	}

	private void readPanel(final String path){
		final LineParser1 p=new LineParser1('\t');
		final boolean[] seen=new boolean[parent.length];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, PANEL_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==5 && p.parseInt(0)==panelGenomes, "Panel order must remain dense and source-bound");
				final int index=byAssembly.get(p.parseString(1));
				require(index>=0 && !seen[index] && genomes.get(index).input.equals(p.parseString(2))
					&& genomes.get(index).canonical.equals(canonical(p.parseString(3))),
					"Every panel genome must bind its actual formal assembly, input ID and FASTA");
				held[root(index)]=true;
				seen[index]=true;
				panelGenomes++;
			}
		}
		require(panelGenomes>0, "An explicit nonempty evaluation panel is required before member voting");
	}

	private void readSources(final String path){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, SOURCE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11, "Current locus sources require their exact eleven-column ledger");
				if(!p.termEquals(member, 1)){continue;}
				final String id=p.parseString(0);
				final int genome=byAssembly.get(p.parseString(2));
				require(genome>=0 && genomes.get(genome).input.equals(p.parseString(3)) && !sourceGenome.contains(id),
					"Every original member locus needs one source-qualified formal assembly");
				sourceGenome.put(id, genome);
			}
		}
		require(sourceGenome.size()>0, "The current source ledger contains no originals for "+member);
	}

	private void readProtected(final String path, final int expected){
		for(final PrepareNcrnaMissingModels.ProtectedPrior prior:PrepareNcrnaMissingModels.readProtectedPriors(path, expected)){
			final int genome=byAssembly.get(prior.assembly);
			require(prior.id.startsWith(member+"_") && genome>=0 && genomes.get(genome).input.equals(prior.input)
				&& genomes.get(genome).canonical.equals(canonical(prior.fasta)) && !sourceGenome.contains(prior.id),
				"Protected source must bind this member, formal assembly/input and canonical original FASTA without a native-ID collision");
			sourceGenome.put(prior.id, genome);
			protectedSources.put(prior.id, prior);
		}
	}

	@Override
	public boolean boundHistoricalSource(final String fullId){
		final PrepareNcrnaMissingModels.ProtectedPrior prior=protectedSources.get(RrnaResourceIO.first(fullId));
		return prior!=null && prior.header.equals(fullId);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Exclusions          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void observe(final String fullId, final byte[] bases){
		require(!filtering, "Observe all original genes before testing any voter for exclusion");
		require(observedSources.add(RrnaResourceIO.first(fullId)),
			"Distinct FASTA headers cannot repeat one original source ID: "+fullId);
		final int genome=genome(fullId);
		final PrepareNcrnaMissingModels.ProtectedPrior prior=protectedSources.get(RrnaResourceIO.first(fullId));
		require(prior==null || prior.header.equals(fullId) && bases.length==prior.length(),
			"Protected FASTA observation must retain its complete historical header and full source length");
		if(!refit){
			if(held[root(genome)]){
				panelSequences.add(bases);
				panelGenes++;
			}else{
				backgroundSequences.add(bases);
				backgroundGenomes.add(genome);
			}
		}
		observedGenes++;
	}

	@Override
	public void finishObservation(){
		require(!filtering && observedGenes>0 && observedGenes==sourceGenome.size()
			&& observedSources.size()==sourceGenome.size(),
			"Every member source-ledger ID must be observed exactly once before any voting or duplicate exclusion");
		for(int gene=0; gene<backgroundSequences.size(); gene++){
			if(panelSequences.contains(backgroundSequences.get(gene))){duplicateGenome[root(backgroundGenomes.get(gene))]=true;}
		}
		backgroundSequences.clear();
		backgroundGenomes.clear();
		filtering=true;
	}

	@Override
	public boolean include(final String fullId, final byte[] bases){
		require(filtering, "All panel sequences must be observed before testing any voter");
		final int genome=genome(fullId);
		if(refit){return true;}
		if(held[root(genome)]){
			excludedAssembly++;
			return false;
		}
		if(panelSequences.contains(bases)){
			duplicateGenome[root(genome)]=true;
			excludedDuplicate++;
			return false;
		}
		return true;
	}

	private int genome(final String fullId){
		final int index=sourceGenome.get(RrnaResourceIO.first(fullId));
		require(index>=0, "Original FASTA IDs must join the current source ledger: "+fullId);
		return index;
	}

	/** Sample only after all voters were tested. Genome-level duplicate exclusion
	 * here is conservative: no padding window from a duplicate-original genome enters
	 * the learned list, while distinct voters elsewhere in that genome remain valid. */
	public void writeSample(final String path, final int count, final long seed){
		require(filtering && panelGenes>0 && count>0 && seed>=0,
			"A training sample follows complete source/sequence exclusions and an explicit sampling seed");
		final int[] eligible=new int[genomes.size()];
		int size=0;
		for(int i=0; i<eligible.length; i++){
			if(!held[root(i)] && !duplicateGenome[root(i)]){eligible[size++]=i;}
		}
		require(size>=count, "Fewer eligible training genomes than requested after whole-assembly and duplicate-voter exclusions: "+size);
		final Random random=new Random(seed);
		for(int i=size-1; i>0; i--){
			final int other=random.nextInt(i+1), value=eligible[i];
			eligible[i]=eligible[other];
			eligible[other]=value;
		}
		RrnaResourceIO.fresh(path, "");
		try(Output out=new Output(path)){
			out.writer.println(PANEL_HEADER);
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<count; i++){
				final Genome g=genomes.get(eligible[i]);
				sampled[eligible[i]]=true;
				out.writer.print(row.clear().append(i).tab().append(g.id).tab().append(g.input)
					.tab().append(g.fasta).tab().append(g.gff).nl());
			}
		}
		sampleCount=count;
		sampleSeed=seed;
	}

	public void writeAudit(final String path){
		require(filtering && (refit || panelGenes>0), "Audit only the completed all-original source join and declared exclusion/refit mode");
		RrnaResourceIO.fresh(path, "");
		try(Output out=new Output(path)){
			out.writer.println("genomes\tpanel_genomes\toriginal_genes\tpanel_alias_genes\tdistinct_panel_sequences\texcluded_assembly_voters\texcluded_duplicate_voters\tunmapped_historical_aliases\ttraining_sample\tsampling_seed");
			out.writer.print(new ByteBuilder().append(genomes.size()).tab().append(panelGenomes).tab().append(observedGenes)
				.tab().append(panelGenes).tab().append(panelSequences.size()).tab().append(excludedAssembly)
				.tab().append(excludedDuplicate).tab().append(unmappedAliases).tab().append(sampleCount).tab().append(sampleSeed).nl());
		}
		RrnaResourceIO.fresh(path+".genomes.tsv", "");
		try(Output out=new Output(path+".genomes.tsv")){
			out.writer.println("assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<genomes.size(); i++){
				final Genome g=genomes.get(i);
				final int group=root(i);
				out.writer.print(row.clear().append(g.id).tab().append(g.input).tab().append(g.canonical)
					.tab().append(g.alias).tab().append(group).tab().append(held[group])
					.tab().append(duplicateGenome[group]).tab().append(sampled[i]).nl());
			}
		}
	}

	public int panelGenomes(){return panelGenomes;}

	/** Read-only original eligibility for other learned statistics after the same
	 * complete observation pass; does not count a voter twice in the profile audit. */
	public boolean allowedOriginal(final String fullId, final byte[] bases){
		require(filtering, "Training statistics require all original sequences observed before exclusion");
		return !held[root(genome(fullId))] && !panelSequences.contains(bases);
	}

	/** Conservative whole-genome capture, excluding duplicate original groups too. */
	public boolean allowedTrainingAssembly(final String assembly){
		require(filtering, "Training genome selection follows complete exact-duplicate observation");
		final int index=byAssembly.get(assembly);
		require(index>=0, "A training assembly must occur in the bound formal manifest: "+assembly);
		final int group=root(index);
		return !held[group] && !duplicateGenome[group];
	}

	private static String canonical(final String path){
		try{return new File(path).getCanonicalPath();}
		catch(IOException e){throw new IllegalArgumentException("Cannot resolve formal assembly identity: "+path, e);}
	}
	private static boolean present(final String text){return text!=null && !text.isEmpty() && !text.equals(".") && !text.equals("NA");}
	private static void require(final boolean condition, final String reason){if(!condition){throw new IllegalArgumentException(reason);}}

	private static final class Genome{
		Genome(final String id_, final String input_, final String fasta_, final String gff_,
				final String accession_, final String canonical_, final int alias_){
			id=id_;
			input=input_;
			fasta=fasta_;
			gff=gff_;
			accession=accession_;
			canonical=canonical_;
			alias=alias_;
		}
		final String id, input, fasta, gff, accession, canonical;
		final int alias;
	}

	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final String member;
	private final boolean refit;
	private final ObjectIntMap<String> byAssembly=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> sourceGenome=new ObjectIntMap<String>(String.class);
	private final ObjectMap<String, PrepareNcrnaMissingModels.ProtectedPrior> protectedSources=
		new ObjectMap<String, PrepareNcrnaMissingModels.ProtectedPrior>(String.class, PrepareNcrnaMissingModels.ProtectedPrior.class);
	private final ObjectSet<String> observedSources=new ObjectSet<String>(String.class);
	private final NcrnaExactSequenceSet panelSequences=new NcrnaExactSequenceSet(true);
	private final ArrayList<byte[]> backgroundSequences=new ArrayList<byte[]>();
	private final IntList backgroundGenomes=new IntList();
	private final int[] parent;
	private final boolean[] held, duplicateGenome, sampled;
	private int panelGenomes, observedGenes, panelGenes, excludedAssembly, excludedDuplicate, unmappedAliases;
	private int sampleCount;
	private long sampleSeed;
	private boolean filtering;
	public static final String MANIFEST_HEADER=NcrnaLearningExclusions.MANIFEST_HEADER;
	public static final String PANEL_HEADER="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
	public static final String SOURCE_HEADER="source_id\tmember\tassembly_id\tinput_id\tcontig\tstart1\tstop1\tstrand\tlength\ttaxid\tchunk";
}
