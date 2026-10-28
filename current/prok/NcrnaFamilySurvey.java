package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import template.ThreadWaiter;

/** Surveys annotated RNA families without calling genes or treating missing annotations as absence.
 * @author Brian Bushnell, Raiden */
public final class NcrnaFamilySurvey {

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFamilySurvey.class, false);
		String map=null, root=null, out=null;int expected=39979, limit=0, threads=8;boolean loci=false;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("map")){map=v;}else if(k.equals("root")){root=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("limit")){limit=Integer.parseInt(v);}
			else if(k.equals("threads") || k.equals("t")){threads=Integer.parseInt(v);}
			else if(k.equals("loci")){loci=Parse.parseBoolean(v);}else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(map!=null && root!=null && out!=null && expected>0 && limit>=0 && threads>0,
			"Specify the complete phylum map, GFF root, expected corpus size and fresh output prefix");
		Shared.setThreads(threads);ReadWrite.setZipThreads(1);ReadWrite.USE_UNPIGZ=false;
		ByteFile.FORCE_MODE_BF1=true;ByteFile.FORCE_MODE_BF2=false;
		final NcrnaFamilySurvey tool=new NcrnaFamilySurvey();
		tool.load(map, root, expected, limit);tool.process(out, threads, loci);
		Shared.closeStream(pp.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Corpus Input         ----------------*/
	/*--------------------------------------------------------------*/

	private void load(final String path, final String root, final int expected, final int limit){
		final ByteFile in=open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final ObjectMap<String, ArrayList<Genome>> pools=new ObjectMap<String, ArrayList<Genome>>(String.class, listClass());
		int total=0;
		try{
			for(byte[] line; (line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==3, "Expected domain/phylum/basename map, no header");
				final String domain=p.parseString(0), phylum=p.parseString(1), name=p.parseString(2);
				require((domain.equals("bacteria") || domain.equals("archaea")) && !phylum.isEmpty()
					&& name.matches("[A-Za-z0-9_.-]+") && names.add(name), "Unique safe corpus names and declared domains required: "+name);
				final String key=domain+":"+phylum;
				ArrayList<Genome> pool=pools.get(key);
				if(pool==null){pool=new ArrayList<Genome>();pools.put(key, pool);}
				pool.add(new Genome(domain, phylum, name, root+"/"+domain+"/"+name+".gff.gz"));total++;
			}
		}finally{require(!in.close(), "Phylum map read failed: "+path);}
		require(total==expected, "Corpus map count changed: "+total+" expected="+expected);
		final ArrayList<String> keys=new ArrayList<String>();
		for(String key:pools.keys()){if(key!=null){keys.add(key);}}
		Collections.sort(keys);
		for(String key:keys){Collections.sort(pools.get(key), (a,b)->a.name.compareTo(b.name));}
		// A bounded canary spans the recorded phyla, rather than taking one taxon's first files.
		for(int rank=0; genomes.size()<total && (limit==0 || genomes.size()<limit); rank++){
			for(String key:keys){
				final ArrayList<Genome> pool=pools.get(key);
				if(rank<pool.size() && (limit==0 || genomes.size()<limit)){genomes.add(pool.get(rank));}
			}
		}
		for(Genome g:genomes){
			final String key=g.domain+":"+g.phylum;int index=phylumIndex.get(key);
			if(index<0){index=phyla.size();phylumIndex.put(key, index);phyla.add(key);}
			g.phylumIndex=index;
		}
		denominators=new int[phyla.size()];results=new GenomeResult[genomes.size()];
		for(Genome g:genomes){denominators[g.phylumIndex]++;}
		require(!genomes.isEmpty(), "Presence fractions require a nonempty selected corpus");
	}

	@SuppressWarnings({"unchecked", "rawtypes"})
	private static Class<ArrayList<Genome>> listClass(){return (Class)ArrayList.class;}

	private void process(final String prefix, final int threads, final boolean detailed){
		RrnaResourceIO.fresh(prefix, ".families.tsv", ".phyla.tsv", ".genomes.tsv", ".loci.tsv", ".summary.tsv", ".aliases.tsv", ".excluded_regulatory.tsv");
		final ByteStreamWriter detail=detailed ? RrnaResourceIO.writer(prefix+".loci.tsv") : null;
		if(detail!=null){detail.println("genome\tfamily\tcontig\tstrand\tintervals1\tlength\tannotation_layers\tfirst_source_line\tlabel");}
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(threads, genomes.size()); i++){workers.add(new Worker(detail));}
		ThreadWaiter.startAndWait(workers);
		if(detail!=null){RrnaResourceIO.finish(detail);}
		for(Worker w:workers){require(w.failure==null && w.success, "GFF worker failed: "+w.failure);}
		final ObjectMap<String, Family> families=new ObjectMap<String, Family>(String.class, Family.class);
		for(Worker w:workers){for(Family f:w.families.values()){
			if(f==null){continue;}Family sum=families.get(f.key);
			if(sum==null){families.put(f.key, f);}else{sum.merge(f);}
		}}
		final ObjectIntMap<String> excluded=new ObjectIntMap<String>(String.class);
		for(Worker w:workers){for(String key:w.excluded.keys()){if(key!=null){addCount(excluded, key, w.excluded.get(key));}}}
		write(prefix, families, excluded);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Per-genome Parsing      ----------------*/
	/*--------------------------------------------------------------*/

	/** Reconstructs RNA branches; a shared gene parent does not merge distinct RNA features. */
	public static ArrayList<Locus> readLoci(final String path, final long[] counters){
		return readLoci(path, counters, new ObjectIntMap<String>(String.class));
	}
	private static ArrayList<Locus> readLoci(final String path, final long[] counters, final ObjectIntMap<String> excluded){
		require(path!=null && counters!=null && counters.length==6, "A GFF scan must retain all six conservation counters");
		final ByteFile in=open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectMap<String, Node> nodes=new ObjectMap<String, Node>(String.class, Node.class);
		final ArrayList<Row> rows=new ArrayList<Row>();long lineNo=0;
		try{
			for(byte[] line; (line=in.nextLine())!=null;){
				lineNo++;
				if(line.length==0){continue;}
				if(line[0]=='#'){
					if(starts(line, "##FASTA")){break;}
					continue;
				}
				p.set(line);require(p.terms()==9, "GFF requires nine fields at "+path+":"+lineNo);counters[0]++;
				if(p.termEquals("CDS", 2)){counters[1]++;continue;}
				final String type=p.parseString(2), attrs=p.parseString(8);
				final String id=attribute(attrs, "ID"), parents=attribute(attrs, "Parent");
				final Node node=node(nodes, id==null ? "@line:"+lineNo : id);
				node.ownPartial|="true".equalsIgnoreCase(attribute(attrs, "partial"))
					|| uncertainRange(attribute(attrs, "start_range")) || uncertainRange(attribute(attrs, "end_range"));
				node.ownPseudo|="true".equalsIgnoreCase(attribute(attrs, "pseudo")) || attribute(attrs, "pseudogene")!=null;
				node.transcript|=type.equals("primary_transcript") || type.equals("transcript");
				if(parents!=null){for(String parent:parents.split(",")){node.parents.add(node(nodes, parent));}}
				final String rfams=rfams(attribute(attrs, "Dbxref"));
				final String gbkey=attribute(attrs, "gbkey"), biotype=attribute(attrs, "gene_biotype");
				node.coding|="CDS".equals(gbkey) || "mRNA".equals(gbkey) || type.equals("mRNA") || "protein_coding".equals(biotype);
				if("CDS".equals(gbkey) || "mRNA".equals(gbkey) || type.equals("mRNA")){continue;}
				final String regulatoryClass=attribute(attrs, "regulatory_class");
				final boolean autonomous=rnaType(type) || rnaType(gbkey) || rnaType(biotype);
				final boolean riboswitch=type.equals("riboswitch") || "riboswitch".equals(regulatoryClass);
				final boolean rna=rfams!=null || autonomous || riboswitch;
				// Generic regulatory annotations include promoters and RBS motifs, not only structured RNAs.
				if(!rna && (regulatoryClass!=null || "regulatory".equals(gbkey) || type.equals("regulatory_region"))){
					addCount(excluded, regulatoryClass==null ? "UNSPECIFIED:"+type : regulatoryClass, 1);continue;
				}
				final boolean carrier=type.equals("gene") || type.equals("primary_transcript") || type.equals("exon");
				if(!rna && !carrier){continue;}
				final int start=p.parseInt(3), stop=p.parseInt(4);
				require(start>=1 && stop>=start && p.parseString(6).length()==1, "Invalid GFF range/strand at "+path+":"+lineNo);
				final String product=attribute(attrs, "product"), note=attribute(attrs, "Note");
				final String label=product!=null ? product : note;
				final int priority=type.equals("exon") ? 3 : type.equals("gene") ? 0 : type.equals("primary_transcript") ? 1 : 2;
				node.primary|=rna && (priority==1 || priority==2);node.geneRna|=rna && priority==0;
				node.orphanRna|=rna;
				node.ownRfams=mergeIds(node.ownRfams, rfams);
				final String effectiveType=carrier && rnaType(gbkey) ? gbkey : carrier && rnaType(biotype) ? biotype : type;
				if(rna && priority!=3){node.ownType=effectiveType;}
				if(label!=null && (node.ownLabel==null || label.compareTo(node.ownLabel)<0)){node.ownLabel=label;}
				rows.add(new Row(node, p.parseString(0), p.parseString(6), start, stop, priority, rna, rfams, label, effectiveType, lineNo));
			}
		}finally{require(!in.close(), "GFF read failed: "+path);}
		for(Node n:nodes.values()){if(n!=null){resolve(n);}}
		for(Node n:nodes.values()){if(n!=null && n.primary){for(Node p0:n.parents){markAncestor(p0);}}}
		final ObjectMap<Node, Component> groups=new ObjectMap<Node, Component>(Node.class, Component.class);
		for(Row row:rows){
			final ArrayList<Node> anchors=anchors(row.node);
			if(anchors.isEmpty()){counters[5]++;continue;}
			for(Node anchor:anchors){
				Component c=groups.get(anchor);
				if(c==null){c=new Component(anchor);groups.put(anchor, c);}c.add(row);
			}
		}
		final ArrayList<Locus> loci=new ArrayList<Locus>();
		final ObjectMap<String, Locus> exact=new ObjectMap<String, Locus>(String.class, Locus.class);
		for(Component c:groups.values()){
			if(c==null || !c.rna){continue;}
			if(c.rfams.size()>0){
				for(String family:c.rfams){makeLocus(c, family, exact, loci, counters);}
			}else{makeLocus(c, c.label==null ? "UNKNOWN:"+c.type : "TEXT:"+c.label, exact, loci, counters);}
		}
		counters[2]=rows.size();counters[3]=loci.size();
		return loci;
	}

	private static void makeLocus(final Component c, final String family, final ObjectMap<String, Locus> exact,
			final ArrayList<Locus> out, final long[] counters){
		final LongList intervals=new LongList();int best=-1, layers=0;String contig=null, strand=null;long first=Long.MAX_VALUE;
		String label=c.label, type=c.type;int labelPriority=-1;
		boolean partial=false, pseudo=false;
		for(Row r:c.rows){
			final String rowFamilies=c.families(r);
			if(c.rfams.size()>1 && (rowFamilies==null || !containsFamily(rowFamilies, family))){continue;}
			if(rowFamilies!=null && !containsFamily(rowFamilies, family)){continue;}
			final int flags=branchFlags(r.node, c.anchor);
			require(flags>=0, "An RNA component row must have a noncoding path to its selected anchor");
			partial|=(flags&1)!=0;pseudo|=(flags&2)!=0;
			if(contig==null){contig=r.contig;strand=r.strand;}
			require(contig.equals(r.contig) && strand.equals(r.strand), "An ID/Parent RNA component crosses contig or strand");
			layers++;first=Math.min(first, r.line);
			final int priority=r.priority==3 ? 1 : r.priority;
			if(r.rna && r.label!=null && (label==null || priority>labelPriority || priority==labelPriority && r.label.compareTo(label)<0)){
				label=r.label;type=r.type;labelPriority=priority;
			}
			if(r.priority<best){continue;}
			if(r.priority>best){intervals.clear();best=r.priority;}
			intervals.add(((long)r.start<<32)|(r.stop&0xffffffffL));
		}
		require(intervals.size>0, "Every family must have an explicit supporting interval: "+family);
		intervals.sort();final LongList merged=new LongList();
		long length=0;int start=(int)(intervals.array[0]>>>32), end=(int)intervals.array[0];
		for(int i=1; i<intervals.size; i++){
			final int s=(int)(intervals.array[i]>>>32), e=(int)intervals.array[i];
			if(s<=(long)end+1){end=Math.max(end, e);}else{merged.add(((long)start<<32)|(end&0xffffffffL));length+=(long)end-start+1;start=s;end=e;}
		}
		merged.add(((long)start<<32)|(end&0xffffffffL));length+=(long)end-start+1;
		require(length>0 && length<=Integer.MAX_VALUE, "RNA union length must fit the exact median length representation");
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<merged.size; i++){if(i>0){b.append(',');}b.append((int)(merged.array[i]>>>32)).append('-').append((int)merged.array[i]);}
		final String spans=b.toString(), signature=family+'\t'+contig+'\t'+strand+'\t'+spans;
		final Locus prior=exact.get(signature);
		if(prior!=null){
			prior.partial|=partial;prior.pseudo|=pseudo;
			prior.earliestSourceLine=Math.min(prior.earliestSourceLine, first);
			prior.sourceIds=mergeIds(prior.sourceIds, c.anchor.id);
			prior.sourceLines=mergeIds(prior.sourceLines, Long.toString(first));counters[4]++;return;
		}
		final Locus locus=new Locus(family, label==null ? "NA" : label, type, contig, strand, spans, (int)length, layers, first);
		locus.partial=partial;locus.pseudo=pseudo;locus.sourceIds=c.anchor.id;locus.sourceLines=Long.toString(first);
		exact.put(signature, locus);out.add(locus);
	}

	private static Node node(final ObjectMap<String, Node> nodes, final String id){
		Node n=nodes.get(id);if(n==null){n=new Node(id);nodes.put(id, n);}return n;
	}
	/** Resolve inherited labels; cycle checks precede branch traversal. */
	private static void resolve(final Node n){
		if(n.state==2){return;}require(n.state==0, "Cyclic ID/Parent ancestry cannot identify RNA loci");n.state=1;
		n.partial=n.ownPartial;n.pseudo=n.ownPseudo;
		String inherited=null, label=null, type=null;
		for(Node p:n.parents){
			resolve(p);if(p.coding){continue;}inherited=mergeIds(inherited, p.resolvedRfams);
			n.partial|=p.partial;n.pseudo|=p.pseudo;
			if(p.resolvedLabel!=null && (label==null || p.resolvedLabel.compareTo(label)<0)){label=p.resolvedLabel;type=p.resolvedType;}
		}
		n.resolvedRfams=n.ownRfams==null ? inherited : n.ownRfams;
		n.resolvedLabel=n.ownLabel==null ? label : n.ownLabel;
		n.resolvedType=n.ownType==null ? type : n.ownType;n.state=2;
	}
	/** A shared exon may reach several RNA anchors. Only paths to the selected
	 * anchor contribute flags; its other RNA parents are siblings, not ancestors. */
	private static int branchFlags(final Node node, final Node anchor){
		assert(node.state==2 && anchor.state==2):"resolve verifies the ancestry DAG before locus reconstruction";
		if(node==anchor){return (anchor.partial ? 1 : 0)|(anchor.pseudo ? 2 : 0);}
		if(node.coding){return -1;}
		int flags=-1;
		for(Node parent:node.parents){
			final int found=branchFlags(parent, anchor);
			if(found>=0){flags=flags<0 ? found : flags|found;}
		}
		return flags<0 ? -1 : flags|(node.ownPartial ? 1 : 0)|(node.ownPseudo ? 2 : 0);
	}
	private static void markAncestor(final Node n){
		if(n.descendantRna){return;}n.descendantRna=true;for(Node p:n.parents){markAncestor(p);}
	}
	/** Keep separate RNA parents, ignore coding parents; common ancestors are metadata, not extra loci. */
	private static ArrayList<Node> anchors(final Node n){
		if(n.anchors!=null){return n.anchors;}n.anchors=new ArrayList<Node>();
		if(n.transcript && n.descendantRna){return n.anchors;}
		if(n.primary){n.anchors.add(n);return n.anchors;}
		if(n.coding || n.geneRna && n.descendantRna){return n.anchors;}
		for(Node p:n.parents){for(Node anchor:anchors(p)){if(!n.anchors.contains(anchor)){n.anchors.add(anchor);}}}
		if(n.anchors.isEmpty() && n.orphanRna && !n.descendantRna){n.anchors.add(n);}
		return n.anchors;
	}
	private static String mergeIds(final String a, final String b){
		if(a==null){return b;}if(b==null || a.equals(b)){return a;}
		final ArrayList<String> ids=new ArrayList<String>();Collections.addAll(ids, a.split(","));
		for(String id:b.split(",")){if(!ids.contains(id)){ids.add(id);}}Collections.sort(ids);return String.join(",", ids);
	}
	private static boolean containsFamily(final String list, final String family){return (","+list+",").contains(","+family+",");}
	/** GFF start/end ranges containing an unknown bound do not establish a full copy. */
	private static boolean uncertainRange(final String range){
		return range!=null && (range.equals(".") || range.startsWith(".,") || range.endsWith(",."));
	}
	private static boolean rnaType(final String type){return type!=null && type.endsWith("RNA") && !type.equals("mRNA");}

	/** GFF attributes split before percent decoding, so escaped delimiters stay in their value. */
	private static String attribute(final String attrs, final String key){
		for(int start=0; start<attrs.length();){
			int stop=attrs.indexOf(';', start);if(stop<0){stop=attrs.length();}
			final int eq=attrs.indexOf('=', start);
			if(eq>start && eq<stop && eq-start==key.length() && attrs.regionMatches(start, key, 0, key.length())){
				return decode(attrs.substring(eq+1, stop));
			}
			start=stop+1;
		}
		return null;
	}
	private static String decode(final String text){
		try{return java.net.URLDecoder.decode(text.replace("+", "%2B"), "UTF-8").replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');}
		catch(java.io.UnsupportedEncodingException e){throw new IllegalStateException(e);}
	}
	private static String rfams(final String db){
		if(db==null){return null;}
		final ArrayList<String> ids=new ArrayList<String>();
		for(String token:db.split(",")){
			if(!token.regionMatches(true, 0, "RFAM:", 0, 5)){continue;}
			final String id=token.substring(5);require(id.matches("RF[0-9]{5}"), "Malformed Rfam accession: "+token);
			if(!ids.contains(id)){ids.add(id);}
		}
		if(ids.isEmpty()){return null;}Collections.sort(ids);return String.join(",", ids);
	}
	private static boolean starts(final byte[] line, final String text){
		if(line.length<text.length()){return false;}for(int i=0; i<text.length(); i++){if(line[i]!=text.charAt(i)){return false;}}return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix, final ObjectMap<String, Family> families, final ObjectIntMap<String> excluded){
		final ByteStreamWriter aliases=RrnaResourceIO.writer(prefix+".aliases.tsv");
		try{
			aliases.println("rfam\tshipped_family");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<SHIPPED_RFAM.length; i++){aliases.print(row.clear().append(SHIPPED_RFAM[i]).tab().append(SHIPPED_FAMILY[i]).nl());}
		}finally{RrnaResourceIO.finish(aliases);}
		final ByteStreamWriter exclusions=RrnaResourceIO.writer(prefix+".excluded_regulatory.tsv");
		try{
			exclusions.println("regulatory_class\tfeatures_excluded_without_Rfam_or_RNA_evidence");
			final ArrayList<String> keys=new ArrayList<String>();for(String key:excluded.keys()){if(key!=null){keys.add(key);}}Collections.sort(keys);
			final ByteBuilder row=new ByteBuilder();
			for(String key:keys){exclusions.print(row.clear().append(key).tab().append(excluded.get(key)).nl());}
		}finally{RrnaResourceIO.finish(exclusions);}
		final ArrayList<Family> sorted=new ArrayList<Family>();
		for(Family f:families.values()){if(f!=null){sorted.add(f);}}
		Collections.sort(sorted, (a,b)->{final int order=Integer.compare(a.block(), b.block());
			if(order!=0){return order;}final int bp=Long.compare(b.bases, a.bases);return bp==0 ? a.key.compareTo(b.key) : bp;});
		final ByteStreamWriter table=RrnaResourceIO.writer(prefix+".families.tsv");
		final ByteStreamWriter byPhylum=RrnaResourceIO.writer(prefix+".phyla.tsv");
		final ByteStreamWriter inputs=RrnaResourceIO.writer(prefix+".genomes.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();long totalLoci=0, rows=0, cds=0, duplicates=0;
		int rfamKeys=0, textKeys=0, unknownKeys=0;
		try{
			table.println("block\tfamily_key\tproduct_or_note\tshipped_family\tterritory_bp\tavg_bp_per_genome\tgenomes_present\tgenomes_total\traw_genome_fraction\tphylum_balanced_genome_fraction\tphyla_present\tphyla_total\tmedian_length\tcopies\tcopies_per_genome\tcopies_per_positive_genome");
			byPhylum.println("family_key\tdomain_phylum\tgenomes_present\tgenomes_total\tpresence_fraction");
			for(Family f:sorted){
				if(f.key.startsWith("RF")){rfamKeys++;}else if(f.key.startsWith("TEXT:")){textKeys++;}else{unknownKeys++;}
				f.lengths.sort();final int n=f.lengths.size;require(n>0 && f.genomes<=genomes.size(), "Counted family must contain valid loci and genome presence");
				final double median=((long)f.lengths.array[(n-1)/2]+f.lengths.array[n/2])/2.0;
				int present=0;double balanced=0;
				for(int i=0; i<phyla.size(); i++){
					require(f.presence[i]<=denominators[i] && denominators[i]>0, "Every phylum uses all sampled genomes, including zero-hit genomes");
					if(f.presence[i]>0){present++;}balanced+=f.presence[i]/(double)denominators[i];
					byPhylum.print(b.clear().append(f.key).tab().append(phyla.get(i)).tab().append(f.presence[i]).tab().append(denominators[i]).tab().append(f.presence[i]/(double)denominators[i], 9).nl());
				}
				balanced/=phyla.size();
				final ArrayList<String> labels=new ArrayList<String>();for(String label:f.labels.keys()){if(label!=null){labels.add(label);}}Collections.sort(labels);
				table.print(b.clear().append(f.block()==0 ? "CANDIDATE" : f.block()==1 ? "SHIPPED" : "UNRESOLVED").tab().append(f.key)
					.tab().append(String.join(" | ", labels)).tab().append(f.shipped).tab().append(f.bases).tab().append(f.bases/(double)genomes.size(), 9)
					.tab().append(f.genomes).tab().append(genomes.size()).tab().append(f.genomes/(double)genomes.size(), 9)
					.tab().append(balanced, 9).tab().append(present).tab().append(phyla.size()).tab().append(median, 1)
					.tab().append(n).tab().append(n/(double)genomes.size(), 9).tab().append(n/(double)f.genomes, 9).nl());
				totalLoci+=n;
			}
			inputs.println("index\tdomain\tphylum\tgenome\tgff\tcompressed_bytes\tmtime_ms\tfeature_rows\tcds_rows\tretained_layers\tcounted_loci\texact_duplicate_components\tlayers_without_noncoding_branch");
			for(int i=0; i<genomes.size(); i++){
				final Genome g=genomes.get(i);final GenomeResult r=results[i];require(r!=null, "Every declared genome must complete, including zero-hit genomes");
				inputs.print(b.clear().append(i).tab().append(g.domain).tab().append(g.phylum).tab().append(g.name).tab().append(g.path).tab().append(r.bytes).tab().append(r.mtime)
					.tab().append(r.counts[0]).tab().append(r.counts[1]).tab().append(r.counts[2]).tab().append(r.counts[3]).tab().append(r.counts[4]).tab().append(r.counts[5]).nl());
				rows+=r.counts[0];cds+=r.counts[1];duplicates+=r.counts[4];
			}
			long genomeLoci=0;for(GenomeResult r:results){genomeLoci+=r.counts[3];}
			require(genomeLoci==totalLoci, "Family copies must conserve the per-genome deduplicated family-locus total");
			require(rfamKeys+textKeys+unknownKeys==sorted.size(), "Rfam, text and unknown keys must partition the annotation-key table");
			summary.println("genomes\tphyla\tannotation_keys\tfeature_rows\tcds_rows\tfamily_loci\texact_duplicate_components\tscope\trfam_keys\ttext_keys\tunknown_keys");
			summary.print(b.clear().append(genomes.size()).tab().append(phyla.size()).tab().append(sorted.size()).tab().append(rows).tab().append(cds)
				.tab().append(totalLoci).tab().append(duplicates).tab().append("ANNOTATION_SURVEY_NOT_BIOLOGICAL_ABSENCE;TERRITORY_SUM_OF_LOCUS_UNIONS;PHYLA_DOMAIN_QUALIFIED")
				.tab().append(rfamKeys).tab().append(textKeys).tab().append(unknownKeys).nl());
		}finally{RrnaResourceIO.finish(table, byPhylum, inputs, summary);}
		System.err.println("NCRNA_SURVEY_PASS genomes="+genomes.size()+" phyla="+phyla.size()+" families="+sorted.size()+" loci="+totalLoci);
	}

	public static String shipped(final String key, final String type, final String label){
		for(int i=0; i<SHIPPED_RFAM.length; i++){if(key.equals(SHIPPED_RFAM[i])){return SHIPPED_FAMILY[i];}}
		if(key.startsWith("RF")){return "NA";}
		if(type.equals("tRNA")){return "tRNA";}if(type.equals("tmRNA")){return "tmRNA";}
		if(type.equals("RNase_P_RNA")){return "RNaseP";}if(type.equals("SRP_RNA")){return "SRP";}
		final String text=label.toLowerCase(Locale.ROOT);
		if(text.equals("16s ribosomal rna") || text.equals("ribosomal rna-16s")){return "16S";}
		if(text.equals("23s ribosomal rna") || text.equals("ribosomal rna-23s")){return "23S";}
		if(text.equals("5s ribosomal rna") || text.equals("ribosomal rna-5s")){return "5S";}if(text.equals("6s rna")){return "6S";}
		return "NA";
	}

	/** Use explicit presence semantics: ObjectIntMap.increment currently adds to an empty-slot sentinel. */
	private static void addCount(final ObjectIntMap<String> map, final String key, final int amount){
		final int old=map.get(key);
		require(amount>0 && (old>=0 || !map.contains(key)), "Exclusion counters contain positive occurrences, never stored negative sentinels: "+key);
		final long sum=(old<0 ? 0L : old)+amount;
		require(sum>0 && sum<=Integer.MAX_VALUE, "Exclusion counts must remain exact without integer overflow: "+key);
		map.put(key, (int)sum);
	}
	private static ByteFile open(final String path){return ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private final class Worker extends Thread{
		Worker(final ByteStreamWriter detail_){detail=detail_;}
		@Override
		public void run(){
			try{
				for(int index; !failed && (index=next.getAndIncrement())<genomes.size();){
					final Genome g=genomes.get(index);final File file=new File(g.path);
					require(file.isFile() && file.length()>0, "Missing GFF: "+g.path);
					final GenomeResult result=new GenomeResult(file.length(), file.lastModified());
					final ArrayList<Locus> loci=readLoci(g.path, result.counts, excluded);
					buffer.clear();
					require(file.length()==result.bytes && file.lastModified()==result.mtime, "GFF changed while reading: "+g.path);
					for(Locus locus:loci){
						Family f=families.get(locus.family);
						if(f==null){f=new Family(locus, phyla.size());families.put(f.key, f);}
						f.add(locus, index, g.phylumIndex);
						if(detail!=null){buffer.append(g.name).tab().append(locus.family).tab().append(locus.contig).tab().append(locus.strand)
							.tab().append(locus.spans).tab().append(locus.length).tab().append(locus.layers).tab().append(locus.line).tab().append(locus.label).nl();}
					}
					if(detail!=null){synchronized(detail){detail.print(buffer);}}
					results[index]=result;
					final int done=completed.incrementAndGet();if(done%1000==0){System.err.println("Genomes completed: "+done+"/"+genomes.size());}
				}
				success=true;
			}catch(Throwable t){failure=t;failed=true;t.printStackTrace(System.err);}
		}
		final ObjectMap<String, Family> families=new ObjectMap<String, Family>(String.class, Family.class);
		final ObjectIntMap<String> excluded=new ObjectIntMap<String>(String.class);
		final ByteStreamWriter detail;final ByteBuilder buffer=new ByteBuilder();
		Throwable failure;boolean success;
	}
	private static final class Genome{
		Genome(final String d, final String p, final String n, final String f){domain=d;phylum=p;name=n;path=f;}
		final String domain, phylum, name, path;int phylumIndex;
	}
	private static final class GenomeResult{
		GenomeResult(final long b, final long m){bytes=b;mtime=m;}
		final long bytes, mtime;final long[] counts=new long[6];
	}
	private static final class Node{
		Node(final String id_){id=id_;}
		final String id;
		boolean ownPartial, ownPseudo, partial, pseudo;
		int state;boolean coding, primary, geneRna, orphanRna, descendantRna, transcript;
		String ownRfams, resolvedRfams, ownLabel, resolvedLabel, ownType, resolvedType;
		final ArrayList<Node> parents=new ArrayList<Node>();
		ArrayList<Node> anchors;
	}
	private static final class Row{
		Row(final Node n, final String c, final String s, final int a, final int z, final int p, final boolean r,
				final String f, final String l, final String t, final long ln){node=n;contig=c;strand=s;start=a;stop=z;priority=p;rna=r;rfams=f;label=l;type=t;line=ln;}
		final Node node;final String contig, strand, label, type;String rfams;final int start, stop, priority;final boolean rna;final long line;
	}
	private static final class Component{
		Component(final Node a){anchor=a;label=a.resolvedLabel;type=a.resolvedType;}
		String families(final Row row){return row.rfams==null ? anchor.resolvedRfams : row.rfams;}
		void add(final Row row){
			rows.add(row);rna|=row.rna;
			final String familyIds=families(row);
			if(familyIds!=null){for(String f:familyIds.split(",")){if(!rfams.contains(f)){rfams.add(f);}}}
			final int metadataPriority=row.priority==3 ? 1 : row.priority;
			if(row.rna && (type==null || row.label!=null && (label==null || metadataPriority>labelPriority
				|| metadataPriority==labelPriority && row.label.compareTo(label)<0))){
				type=row.type;label=row.label;labelPriority=metadataPriority;
			}
		}
		final ArrayList<Row> rows=new ArrayList<Row>();final ArrayList<String> rfams=new ArrayList<String>();
		final Node anchor;
		boolean rna;String label, type;int labelPriority=-1;
	}
	public static final class Locus{
		Locus(final String f, final String l, final String t, final String c, final String s, final String spans_, final int n, final int layers_, final long ln){
			family=f;label=l;type=t;contig=c;strand=s;spans=spans_;length=n;layers=layers_;line=ln;earliestSourceLine=ln;
		}
		public final String family, label, type, contig, strand, spans;public final int length, layers;public final long line;
		/** Source flags are inherited through noncoding parents and unioned across exact duplicates. */
		public boolean partial, pseudo;
		public String sourceIds, sourceLines;
		/** Deterministic provenance across exact duplicate components; legacy line is unchanged. */
		public long earliestSourceLine;
	}
	private static final class Family{
		Family(final Locus l, final int phyla){key=l.family;shipped=shipped(l.family, l.type, l.label);presence=new int[phyla];}
		void add(final Locus l, final int genome, final int phylum){
			assert(key.equals(l.family)):"A family accumulator cannot mix Rfam accessions";
			labels.add(l.label);bases+=l.length;lengths.add(l.length);
			if(lastGenome!=genome){lastGenome=genome;genomes++;presence[phylum]++;}
		}
		void merge(final Family f){
			require(key.equals(f.key) && shipped.equals(f.shipped), "Family metadata must agree across workers: "+key);
			bases+=f.bases;genomes+=f.genomes;lengths.addAll(f.lengths);
			for(String label:f.labels.keys()){if(label!=null){labels.add(label);}}for(int i=0; i<presence.length; i++){presence[i]+=f.presence[i];}
		}
		int block(){return !shipped.equals("NA") ? 1 : key.startsWith("RF") ? 0 : 2;}
		final String key, shipped;final ObjectSet<String> labels=new ObjectSet<String>(String.class);
		final IntList lengths=new IntList();final int[] presence;long bases;int genomes, lastGenome=-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<String> phyla=new ArrayList<String>();
	private final ObjectIntMap<String> phylumIndex=new ObjectIntMap<String>(String.class);
	private final AtomicInteger next=new AtomicInteger(), completed=new AtomicInteger();
	private int[] denominators;private GenomeResult[] results;
	private volatile boolean failed=false;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Shipped families, including eukaryotic exclusions. RF02542 is nuclear SSU, not 16S.
	private static final String[] SHIPPED_RFAM={"RF00177", "RF01959", "RF02540", "RF02541", "RF00001", "RF00005",
		"RF00010", "RF00011", "RF00373", "RF00023", "RF01849", "RF01850", "RF01851", "RF04321", "RF00013", "RF01685",
		"RF00169", "RF04183", "RF01854", "RF01857", "RF01852", "RF01960", "RF02542", "RF00002", "RF02543"};
	private static final String[] SHIPPED_FAMILY={"16S", "16S", "23S", "23S", "5S", "tRNA", "RNaseP", "RNaseP", "RNaseP",
		"tmRNA", "tmRNA", "tmRNA", "tmRNA", "tmRNA", "6S", "6S", "SRP", "SRP", "SRP", "SRP", "tRNA", "18S", "18S", "5.8S", "eukLSU"};
}
