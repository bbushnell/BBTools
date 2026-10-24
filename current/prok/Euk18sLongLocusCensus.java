package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.function.IntFunction;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import tax.TaxNode;
import tax.TaxTree;

/** Count long CM-included 18S loci in saved truth ledgers; never reads sequences.
 * Length is the genomic hit span, not mature RNA length or an intron diagnosis.
 * @author Ganyu
 */
public final class Euk18sLongLocusCensus {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,Euk18sLongLocusCensus.class,false);
		String manifest=null,treePath=null,prefix=null;
		for(String arg:pp.args){
			final int split=arg.indexOf('=');require(split>0,"Expected flag=value");
			final String key=arg.substring(0,split).toLowerCase(Locale.ROOT),value=arg.substring(split+1);
			if(key.equals("in")){manifest=value;}else if(key.equals("tree")){treePath=value;}
			else if(key.equals("outprefix")){prefix=value;}else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		require(manifest!=null && treePath!=null && prefix!=null,"Declare saved-ledger manifest, taxonomy and fresh outprefix");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;
		final TaxTree tree=TaxTree.loadTaxTreeFile(treePath);
		new Euk18sLongLocusCensus().run(manifest,prefix,id->tree.getNode(id,true));
		Shared.closeStream(pp.outstream);
	}

	public void run(String manifest,String prefix,IntFunction<TaxNode> lookup){
		require(lookup!=null,"Taxonomy lookup supplies organism labels; it is not inferred from filename text");
		RrnaResourceIO.fresh(prefix,".summary.tsv",".organisms.tsv",".long_loci.tsv.gz");
		final ByteFile input=RrnaResourceIO.open(manifest);
		final ByteStreamWriter longLoci=RrnaResourceIO.writer(prefix+".long_loci.tsv.gz");
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> clades=new ObjectSet<String>(String.class);
		longLoci.println(RrnaSpeciesKeys.TRUTH_HEADER);
		try{
			RrnaResourceIO.header(input,"clade\ttruth\tpopulation_summary",manifest);
			for(byte[] row;(row=input.nextLine())!=null;){
				p.set(row);require(p.terms()==3,"Manifest requires clade/truth/population_summary");
				final String clade=p.parseString(0);
				require(validClade(clade) && clades.add(clade),"Each recognized clade may appear only once: "+clade);
				final Counts counts=new Counts(clade);
				readTruth(clade,p.parseString(1),counts,longLoci,lookup);
				checkPopulation(p.parseString(2),counts);
				summaries.add(counts);pooled.add(counts);
			}
		}finally{final boolean error=input.close();RrnaResourceIO.finish(longLoci);require(!error,"Input manifest read failed");}
		require(clades.size()>0,"An empty manifest is not a corpus census");
		summaries.add(pooled);summaries.sort((a,b)->a.clade.compareTo(b.clade));
		groups.sort((a,b)->a.key.compareTo(b.key));
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv"),organisms=RrnaResourceIO.writer(prefix+".organisms.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			summary.println("clade\tpopulation\tloci\tgt2200\tgt2500\tmedianGenomicSpan\tp99GenomicSpan\tp998GenomicSpan\tmaxGenomicSpan");
			organisms.println("clade\tassembly\tsource_taxid\ttaxon_name\tspecies_taxid\tspecies_name\tpopulation\tloci\tgt2200\tgt2500");
			for(Counts c:summaries){for(int population=0;population<3;population++){
				appendCounts(b.clear().append(c.clade).tab(),c,population,true);summary.print(b);
			}}
			for(Group g:groups){
				if(g.counts.values[0][1]==0){continue;}
				final TaxNode taxon=lookup.apply(g.taxid);require(taxon!=null,"Missing source taxid in pinned taxonomy: "+g.taxid);
				final int speciesId=RrnaSpeciesKeys.species(g.taxid,lookup);
				final TaxNode species=speciesId>0?lookup.apply(speciesId):null;
				for(int population=0;population<3;population++){
					b.clear().append(g.clade).tab().append(g.assembly).tab().append(g.taxid).tab().append(label(taxon)).tab()
						.append(speciesId).tab().append(species==null?"UNRESOLVED":label(species)).tab();
					appendCounts(b,g.counts,population,false);organisms.print(b);
				}
			}
		}finally{RrnaResourceIO.finish(summary,organisms);}
		System.err.println("LONG_LOCUS_CENSUS_PASS clades="+clades.size()+" included="+pooled.values[0][0]
			+" gt2200="+pooled.values[0][1]+" gt2500="+pooled.values[0][2]+" sequenceFilesOpened=0 alignments=0");
	}

	void readTruth(String clade,String path,Counts counts,ByteStreamWriter longLoci,IntFunction<TaxNode> lookup){
		final ByteFile input=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		Group previous=null;
		try{
			RrnaResourceIO.header(input,RrnaSpeciesKeys.TRUTH_HEADER,path);
			for(byte[] row;(row=input.nextLine())!=null;){
				p.set(row);require(p.terms()==22 && p.termEquals(clade,1),"Truth schema/clade differs from manifest: "+path);
				final int taxid=p.parseInt(3),start=p.parseInt(5),stop=p.parseInt(6),length=p.parseInt(8),contigLength=p.parseInt(9);
				require(taxid>0 && start>0 && stop>=start && stop<=contigLength && length==(long)stop-start+1
					&& (p.termEquals("+",7)||p.termEquals("-",7)),"Invalid inclusive truth span/taxid: "+path);
				require(seen.add(p.parseString(0)),"Duplicate locus ID would double-count corpus copies: "+path);
				final boolean admitted=flag(p,17),construction=flag(p,21);
				require(construction==(admitted && p.termEquals("no_keyword",20)),"Construction population must retain the declared admitted nuclear-only rule");
				Euk18sM1Grade.validComp(p.parseString(20));
				if(previous==null || previous.taxid!=taxid || !p.termEquals(previous.assembly,2)){
					final String assembly=p.parseString(2),key=clade+'\t'+assembly+'\t'+taxid;
					previous=byAssembly.get(key);
					if(previous==null){previous=new Group(clade,assembly,taxid,key);byAssembly.put(key,previous);groups.add(previous);}
				}
				counts.observe(length,admitted,construction);previous.counts.observe(length,admitted,construction);
				if(length>2200){longLoci.println(row);}
			}
		}finally{require(!input.close(),"Truth input read failed: "+path);}
		require(counts.values[0][0]>0,"Declared clade has no included loci: "+clade);
	}

	static void checkPopulation(String path,Counts observed){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final boolean[] seen=new boolean[3];
		try{
			RrnaResourceIO.header(in,"population\tloci\tmin\tp01\tp05\tp50\tp95\tp99\tmax",path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==9,"Population summary schema changed: "+path);
				int index=-1;for(int i=0;i<3;i++){if(p.termEquals(POPULATIONS[i],0)){index=i;}}
				require(index>=0 && !seen[index],"Unknown/repeated population: "+path);seen[index]=true;
				require(p.parseLong(1)==observed.values[index][0],"Truth count differs from sealed population summary: "+path);
			}
		}finally{require(!in.close(),"Population summary read failed");}
		for(boolean value:seen){require(value,"Missing included/admitted/construction summary");}
	}
	static boolean flag(LineParser1 p,int column){require(p.termEquals("0",column)||p.termEquals("1",column),"Population flag must be0 or1");return p.termEquals("1",column);}
	static boolean validClade(String c){for(String value:CLADES){if(value.equals(c)){return true;}}return false;}
	static String label(TaxNode node){require(node.name!=null && node.name.indexOf('\t')<0 && node.name.indexOf('\n')<0,"Taxonomy name must fit a TSV field");return node.name;}
	static void appendCounts(ByteBuilder b,Counts c,int population,boolean quantiles){
		final long[] n=c.values[population];require(n[0]>=n[1] && n[1]>=n[2],"Nested strict length thresholds must conserve counts");
		b.append(POPULATIONS[population]).tab().append(n[0]).tab().append(n[1]).tab().append(n[2]);
		if(quantiles){
			require(c.histograms!=null,"Genomic-span quantiles require the multiplicity-preserving histogram");
			final long[] h=c.histograms[population];
			b.tab().append(quantile(h,n[0],1,2)).tab().append(quantile(h,n[0],99,100)).tab()
				.append(quantile(h,n[0],998,1000)).tab().append(quantile(h,n[0],1,1));
		}
		b.nl();
	}
	public static int quantile(long[] histogram,long count,long numerator,long denominator){
		require(count>=0 && numerator>0 && numerator<=denominator,"Nearest rank requires a valid fraction and observed count");
		if(count==0){return -1;}
		final long rank=Math.addExact(Math.multiplyExact(count,numerator),denominator-1)/denominator;
		long sum=0;for(int i=1;i<histogram.length;i++){sum+=histogram[i];if(sum>=rank){return i;}}
		throw new IllegalStateException("Length histogram does not conserve observed loci");
	}
	static void require(boolean value,String why){if(!value){throw new IllegalArgumentException(why);}}
	public static final class Counts {
		public Counts(String clade){this(clade,true);}
		Counts(String clade,boolean histogram){this.clade=clade;histograms=histogram?new long[3][4096]:null;}
		public void observe(int length,boolean admitted,boolean construction){
			require(length>0 && (!construction||admitted),"Construction must be a subset of admitted positive-length loci");
			for(int i=0;i<3;i++){
				if(i==1&&!admitted || i==2&&!construction){continue;}
				values[i][0]++;if(length>2200){values[i][1]++;}if(length>2500){values[i][2]++;}
				if(histograms!=null){if(length>=histograms[i].length){histograms[i]=Arrays.copyOf(histograms[i],length+1);}histograms[i][length]++;}
			}
		}
		public void add(Counts other){
			require(histograms!=null && other.histograms!=null,"Pooled quantiles must merge clade histograms, not average quantiles");
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){values[i][j]+=other.values[i][j];}
				if(histograms[i].length<other.histograms[i].length){histograms[i]=Arrays.copyOf(histograms[i],other.histograms[i].length);}
				for(int j=0;j<other.histograms[i].length;j++){histograms[i][j]+=other.histograms[i][j];}
			}
		}
		final String clade;public final long[][] values=new long[3][3],histograms;
	}
	static final class Group {
		Group(String clade,String assembly,int taxid,String key){this.clade=clade;this.assembly=assembly;this.taxid=taxid;this.key=key;counts=new Counts(clade,false);}
		final String clade,assembly,key;final int taxid;final Counts counts;
	}
	static final String[] POPULATIONS={"included","admitted","construction"};
	public static final String[] CLADES={"fungi","protozoa","plant","invertebrate","vertebrate_mammalian","vertebrate_other"};
	final Counts pooled=new Counts("pooled_euk");
	final ArrayList<Counts> summaries=new ArrayList<Counts>();
	final ObjectMap<String,Group> byAssembly=new ObjectMap<String,Group>(String.class,Group.class);
	final ArrayList<Group> groups=new ArrayList<Group>();
}
