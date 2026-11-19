package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.LongIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaReferenceEvidence.Evidence;
import prok.NcrnaReferenceEvidence.Hit;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import structures.ByteBuilder;

/** Imports one completed batch at a time, retaining real native report witnesses.
 * Preparation fingerprints are declared execution provenance, not hashes recomputed
 * by this reader. Actual native command footers bind model, input and original-GA
 * inclusion. Complete captured maps replace redundant current-FASTA observations.
 * Native completion counts are batch/model counts. Per-input counts are derived
 * only after all native model reports have been validated, including empty reports.
 * No per-input native search, CPU measurement or native-line renumbering is invented.
 * One worker owns this mutable reader and its bounded batch of Evidence objects.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceBatch{

	/** Loads immutable preparation membership before workers acquire sequence evidence. */
	static Task task(final String id, final String prepared, final String results,
		final long expectedHits, final NcrnaEvidenceCatalog catalog){
		return task(id, prepared, results, results, expectedHits, catalog);
	}

	/** Opens recovered bytes separately from the immutable original native path. */
	static Task task(final String id, final String prepared, final String results, final String physical,
		final long expectedHits, final NcrnaEvidenceCatalog catalog){
		return new Task(id, prepared, results, physical, expectedHits, catalog);
	}

	/** The bundled collector has scan-local batch IDs and no fabricated hit total. */
	static Task bundledTask(final String group, final String id, final String prepared, final String results,
		final String bundle, final long genomes, final long bases, final long records, final NcrnaEvidenceCatalog catalog){
		final Task task=new Task(id, prepared, results, results, -1, catalog, group, bundle);
		require(task.genomes.size()==genomes && task.bases==bases && task.records==records,
			"Bundled collector dimensions must match the native input receipt and preparation");
		return task;
	}

	/** Check each CM header once before workers start; HMM-filter headers are separate.
	 * This is metadata binding, not a substitute for the recipe's byte authentication.
	 */
	static void validateBundle(final String path, final Roster roster){
		final Words words=new Words();
		final boolean[] seen=new boolean[roster.queried()];
		boolean inModel=false;
		String name=null, accession=null;
		long length=-1;
		double ga=Double.NaN;
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				words.set(row);
				if(words.eq(0, "INFERNAL1/a")){
					require(!inModel, "Unterminated CM in bundle");
					inModel=true;
					name=null;
					accession=null;
					length=-1;
					ga=Double.NaN;
				}else if(inModel && words.eq(0, "NAME")){
					require(name==null, "Repeated CM name");
					name=words.text(1);
				}else if(inModel && words.eq(0, "ACC")){
					require(accession==null, "Repeated CM accession");
					accession=words.text(1);
				}else if(inModel && words.eq(0, "CLEN")){
					require(length<0, "Repeated CM length");
					length=words.number(1, 0, 0);
				}else if(inModel && words.eq(0, "GA")){
					require(Double.isNaN(ga), "Repeated CM GA threshold");
					ga=Double.parseDouble(words.text(1));
				}else if(inModel && words.eq(0, "//")){
					final Model model=accession==null ? null : roster.model(accession);
					require(model!=null && model.queryIndex()>=0 && !seen[model.queryIndex()]
						&& model.queryName().equals(name) && model.length==length && model.ga==ga,
						"Bundled CM metadata differs from the captured queried roster: "+accession);
					seen[model.queryIndex()]=true;
					inModel=false;
				}
			}
		}
		require(!inModel, "Unterminated final CM in bundle");
		for(final boolean present:seen){require(present, "Every queried model must occur in the bound bundle");}
	}

	/** Reads each model table once; successful return conserves the whole batch. */
	static Evidence[] read(final Task task, final NcrnaEvidenceCatalog catalog,
		final Roster roster, final double searchSpaceMb){
		final NcrnaReferenceBatch reader=new NcrnaReferenceBatch(task, catalog, roster, searchSpaceMb);
		reader.readMap();
		reader.readCompletion();
		long total=0;
		if(task.bundle!=null){
			reader.readBundledSearch(new File(task.physicalResults, "search.out").getPath());
			total=reader.readHits(new File(task.results, "hits.tbl").getPath(),
				new File(task.physicalResults, "hits.tbl").getPath(), null, Double.NaN);
			for(int i=0; i<roster.queried(); i++){
				require(reader.observed[i]==reader.expected[i], "Bundled native hit count differs for "+roster.query(i).accession);
			}
		}else{
			for(int i=0; i<roster.queried(); i++){
				final Model model=roster.query(i);
				final String point=new File(task.results, model.accession).getPath();
				final String physical=new File(task.physicalResults, model.accession).getPath();
				require(line(new File(physical, "exit.txt").getPath()).equals("0"), "A failed model search cannot supply empty evidence: "+point);
				reader.readSearch(new File(physical, "search.out").getPath(), model, reader.expected[i]);
				final long observed=reader.readHits(new File(point, "hits.tbl").getPath(), new File(physical, "hits.tbl").getPath(), model, reader.reportBits[i]);
				require(observed==reader.expected[i], "Native reported-hit count differs from the complete table: "+model.accession);
				total=Math.addExact(total, observed);
			}
		}
		require(task.bundle!=null || total==task.expectedHits, "Batch collector and actual native hit totals differ");
		long distributed=0;
		for(final Evidence evidence:reader.evidence){distributed=Math.addExact(distributed, evidence.size());}
		require(distributed==total, "Every native hit must map to exactly one input");
		return reader.evidence;
	}

	private NcrnaReferenceBatch(final Task task_, final NcrnaEvidenceCatalog catalog,
		final Roster roster_, final double searchSpaceMb_){
		task=task_;
		roster=roster_;
		searchSpaceMb=searchSpaceMb_;
		require(Double.isFinite(searchSpaceMb) && searchSpaceMb>0, "Native -Z binding requires a positive finite search space");
		evidence=new Evidence[task.genomes.size()];
		for(int i=0; i<evidence.length; i++){
			final Genome genome=task.genomes.get(i);
			final Observation captured=catalog.captured(genome, task.mapping, task.mapSha80);
			evidence[i]=new Evidence(null, captured, roster);
			inputIndex.put(Long.parseLong(genome.id), i);
		}
		expected=new long[roster.queried()];
		observed=new long[roster.queried()];
		reportBits=new double[roster.queried()];
	}

	/*--------------------------------------------------------------*/
	/*----------------       Captured Targets       ----------------*/
	/*--------------------------------------------------------------*/

	private void readMap(){
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder name=new ByteBuilder();
		long count=0;
		try(Input in=new Input(task.mapping)){
			RrnaResourceIO.header(in.reader, PrepareCmSearchBatches.TARGET_HEADER, task.mapping);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==6 && p.termEquals(task.id, 0), "Captured map contains a foreign batch or invalid schema");
				final int input=inputIndex.get(p.parseLong(2));
				require(input>=0 && p.termEquals(evidence[input].observed.genome.id, 2), "Captured target names an input outside this batch");
				final Observation observation=evidence[input].observed;
				final long ordinal=p.parseLong(3);
				final int length=p.parseInt(5);
				name.clear().append('i').append(observation.genome.id).append('c').append(ordinal);
				require(same(p, 1, name), "Synthetic target does not encode its exact input and contig ordinal");
				probe.set(p, 1);
				require(!targets.contains(probe), "Repeated synthetic target in captured map");
				final String original=PrepareCmSearchBatches.targetToken(text(p, 4));
				require(p.termEquals(original, 4), "Captured original target must be one complete literal token");
				observation.addCaptured(original, length, ordinal);
				targets.put(probe.copy(), new Mapping(input, p.parseString(1), original));
				count++;
			}
		}
		require(count==task.records, "Captured target map does not cover all prepared contigs");
		for(final Evidence item:evidence){item.observed.finishCaptured();}
	}

	/** Reads actual batch/model completion, not inferred per-genome native receipts. */
	private void readCompletion(){
		require(line(new File(task.physicalResults, "exit.txt").getPath()).equals("0"), "Incomplete batch search cannot be imported");
		final boolean bundled=task.bundle!=null;
		if(!bundled){
			final byte[] copied=oneRow(new File(task.physicalResults, "input.tsv").getPath(), BINDING_HEADER);
			require(Arrays.equals(copied, task.binding), "Search must retain the exact accepted preparation binding");
		}
		final LineParser1 p=new LineParser1('\t');
		final boolean[] seen=new boolean[roster.queried()];
		int count=0;
		final String path=new File(task.physicalResults, "completion.tsv").getPath();
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, bundled ? BUNDLED_COMPLETION_HEADER : COMPLETION_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==(bundled ? 7 : 9) && p.termEquals(task.id, 0)
					&& p.termEquals("COMPLETE", bundled ? 6 : 8), "Native completion must identify this complete batch");
				int query=-1;
				for(int i=0; i<roster.queried(); i++){if(p.termEquals(roster.query(i).accession, 1)){query=i;}}
				require(query>=0 && !seen[query], "Unknown or repeated completed model");
				require(p.parseLong(2)==task.records && p.parseLong(3)==task.bases && p.parseLong(4)==Math.multiplyExact(task.bases, 2),
					"Native batch dimensions differ from the exact prepared population");
				expected[query]=p.parseLong(5);
				require(expected[query]>=0, "Native reported-hit count cannot be negative");
				if(!bundled){
					reportBits[query]=p.parseDouble(7);
					require(p.parseDouble(6)==roster.query(query).ga && Double.isFinite(reportBits[query])
						&& reportBits[query]>0 && reportBits[query]<roster.query(query).ga, "Invalid native hit count or original-GA/report binding");
				}
				seen[query]=true;
				count++;
			}
		}
		require(count==seen.length, "Missing model completion cannot become zero-hit evidence");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Native Evidence        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads physical bytes but keeps original path/line witnesses and footer checks. */
	private long readHits(final String path, final String physicalPath, final Model selected, final double lower){
		final RrnaCmTruthBuilder.Fields fields=new RrnaCmTruthBuilder.Fields();
		long number=0, count=0;
		boolean program=false, version=false, command=false, finished=false;
		Model previous=selected;
		try(Input in=new Input(physicalPath)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				number++;
				if(row.length==0){continue;}
				if(finished){throw new IllegalArgumentException("Content follows native completion: "+path+':'+number);}
				if(row[0]=='#'){
					final String text=new String(row, StandardCharsets.UTF_8);
					if(text.startsWith("# Program:")){
						require(!program && text.substring(10).trim().equals("cmsearch"), "Expected one native cmsearch program marker");
						program=true;
					}else if(text.startsWith("# Version:")){
						require(!version && text.substring(10).trim().startsWith("1.1.5 "), "Expected recorded Infernal 1.1.5 format");
						version=true;
					}else if(text.startsWith("# Option settings:")){
						require(!command, "Repeated native command footer");
						checkCommand(text.substring(18).trim(), selected, lower, path);
						command=true;
					}else if(text.equals("# [ok]")){finished=true;}
					continue;
				}
				fields.set(row);
				if(selected==null && (previous==null || !fields.eq(3, previous.accession))){previous=roster.model(fields.text(3));}
				final Model model=previous;
				if(model==null || model.queryIndex()<0 || !fields.eq(4, "cm") || !fields.eq(3, model.accession) || !fields.eq(2, model.queryName())){
					throw new IllegalArgumentException("Unexpected model/type in native report: "+path+':'+number);
				}
				probe.set(row, fields.start[0], fields.end[0]);
				final Mapping target=targets.get(probe);
				if(target==null){throw new IllegalArgumentException("Unknown synthetic native target: "+path+':'+number);}
				final Evidence item=evidence[target.input];
				final Hit hit=new Hit(path, number, row, target.original, target.synthetic, model, fields, item.observed, true);
				require(selected!=null || hit.included, "A --cut_ga bundled report cannot contain a natively excluded hit");
				if(!hit.range.equals("VALID_INTERVAL")){throw new IllegalArgumentException("Native coordinates leave the captured contig: "+path+':'+number);}
				item.addBatchHit(hit);
				observed[model.queryIndex()]++;
				count++;
			}
		}
		require(program && version && command && finished, "Native report lacks program/version/bound-command/completion footer: "+path);
		return count;
	}

	/** The real native footer, rather than a manually reconstructed command, binds inclusion. */
	private void checkCommand(final String command, final Model model, final double lower, final String table){
		if(model==null){CompareCmSearchCosts.checkCommand(command, task.bundle, Double.NaN, Double.NaN, false);}
		else{CompareCmSearchCosts.checkCommand(command, model.path, model.ga, lower, true);}
		final String[] argv=command.trim().split(" +");
		require(argv[argv.length-1].equals(task.fasta), "Native command searched a different prepared batch");
		int cpu=0, space=0, format=0, output=0, searchReport=0, noHmm=0, noAli=0;
		for(int i=0; i<argv.length-2; i++){
			final String option=argv[i];
			if(option.equals("--cpu")){
				require(i+1<argv.length-2 && argv[++i].equals("1"), "Expected one native worker per model search");
				cpu++;
			}else if(option.equals("-Z")){
				require(i+1<argv.length-2 && Double.parseDouble(argv[++i])==searchSpaceMb, "Native full search space differs from import binding");
				space++;
			}else if(option.equals("--tformat")){
				require(i+1<argv.length-2 && argv[++i].equals("fasta"), "Native input format must be FASTA");
				format++;
			}else if(option.equals("--tblout")){
				require(i+1<argv.length-2 && argv[++i].equals(table), "Native command names a different hit table");
				output++;
			}else if(option.equals("-o")){
				final String report=new File(new File(table).getParentFile(), "search.out").getPath();
				require(i+1<argv.length-2 && argv[++i].equals(report), "Native command names a different search report");
				searchReport++;
			}else if(option.equals("--nohmmonly")){noHmm++;}
			else if(option.equals("--noali")){noAli++;}
		}
		require(cpu==1 && space==1 && format==1 && output==1 && searchReport==1 && noHmm==1 && noAli==1, "Native command is missing or repeats a required search option");
	}

	/** Independently checks the native search report's query, dimensions and hit count. */
	private void readSearch(final String path, final Model model, final long hits){
		final Words words=new Words();
		int queries=0, accessions=0, dimensions=0, totals=0, ends=0;
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				words.set(row);
				require(ends==0 || words.count==0, "Content follows native search completion");
				if(words.eq(0, "Query:")){
					require(words.eq(1, model.queryName()) && words.clen(2)==model.length, "Native query name/length differs from the model roster");
					queries++;
				}else if(words.eq(0, "Accession:")){
					require(words.eq(1, model.accession), "Native query accession differs from the model roster");
					accessions++;
				}else if(words.eq(0, "Target") && words.eq(1, "sequences:")){
					require(words.count>3 && row[words.start[3]]=='(', "Malformed native searched-residue count");
					require(words.number(2, 0, 0)==task.records && words.number(3, 1, 0)==Math.multiplyExact(task.bases, 2), "Native report searched different target dimensions");
					dimensions++;
				}else if(words.eq(0, "Total") && words.eq(1, "CM") && words.eq(2, "hits") && words.eq(3, "reported:")){
					require(words.number(4, 0, 0)==hits, "Native report and completion hit counts differ");
					totals++;
				}else if(words.eq(0, "[ok]")){ends++;}
			}
		}
		require(queries==1 && accessions==1 && dimensions==1 && totals==1 && ends==1, "Native search report is incomplete or contains multiple queries");
	}

	/** Each query, including a zero-hit query, must have its own complete accounting. */
	private void readBundledSearch(final String path){
		final Words words=new Words();
		final boolean[] seen=new boolean[roster.queried()];
		int query=-1;
		int accessions=0, dimensions=0, totals=0, ends=0;
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				words.set(row);
				require(ends==0 || words.count==0, "Content follows native search completion");
				if(words.eq(0, "Query:")){
					require(query<0 || accessions==1 && dimensions==1 && totals==1, "Incomplete bundled query accounting");
					query=-1;
					for(int i=0; i<roster.queried(); i++){
						if(words.eq(1, roster.query(i).queryName())){
							require(query<0, "Query name must identify one roster model");
							query=i;
						}
					}
					require(query>=0 && !seen[query] && words.clen(2)==roster.query(query).length,
						"Unknown, repeated or length-mismatched bundled query");
					seen[query]=true;
					accessions=0;
					dimensions=0;
					totals=0;
				}else if(words.eq(0, "Accession:")){
					require(query>=0 && words.eq(1, roster.query(query).accession), "Bundled query accession differs from roster");
					accessions++;
				}else if(words.eq(0, "Target") && words.eq(1, "sequences:")){
					require(query>=0 && words.count>3 && row[words.start[3]]=='('
						&& words.number(2, 0, 0)==task.records && words.number(3, 1, 0)==Math.multiplyExact(task.bases, 2),
						"Bundled native query searched different target dimensions");
					dimensions++;
				}else if(words.eq(0, "Total") && words.eq(1, "CM") && words.eq(2, "hits") && words.eq(3, "reported:")){
					require(query>=0 && words.number(4, 0, 0)==expected[query], "Bundled report and completion hit counts differ");
					totals++;
				}else if(words.eq(0, "[ok]")){ends++;}
			}
		}
		require(query>=0 && accessions==1 && dimensions==1 && totals==1 && ends==1, "Bundled native search report is incomplete");
		for(final boolean present:seen){require(present, "Missing native query cannot be treated as zero-hit coverage");}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Immutable Binding      ----------------*/
	/*--------------------------------------------------------------*/

	static final class Task{
		Task(final String id_, final String prepared_, final String results_, final String physicalResults_,
			final long expectedHits_, final NcrnaEvidenceCatalog catalog){
			this(id_, prepared_, results_, physicalResults_, expectedHits_, catalog, null, null);
		}
		Task(final String id_, final String prepared_, final String results_, final String physicalResults_,
			final long expectedHits_, final NcrnaEvidenceCatalog catalog, final String group_, final String bundle_){
			id=id_;
			bundle=bundle_;
			key=bundle==null ? id : group_+":"+id;
			prepared=prepared_;
			results=results_;
			physicalResults=physicalResults_;
			expectedHits=expectedHits_;
			require(id.matches("0|[1-9][0-9]*") && new File(prepared).isAbsolute() && new File(results).isAbsolute()
				&& new File(physicalResults).isAbsolute() && (bundle==null ? expectedHits>=0
					: expectedHits==-1 && group_!=null && group_.matches("[A-Za-z0-9_.-]+") && new File(bundle).isAbsolute()),
				"Batch collector needs a canonical ID, absolute directories and nonnegative hit count");
			final LineParser1 p=new LineParser1('\t');
			final long count;
			if(bundle!=null){
				binding=null;
				fasta=new File(prepared, "batch.fa").getPath();
				mapping=new File(prepared, "targets.tsv").getPath();
				members=new File(prepared, "members.tsv").getPath();
				final String[] hashes=readInputReceipt(new File(physicalResults, "input.tsv").getPath(), prepared);
				fastaSha80=hashes[0];
				mapSha80=hashes[1];
				membersSha80=hashes[2];
				p.set(oneRow(new File(prepared, "summary.tsv").getPath(), PREPARED_HEADER));
				require(p.terms()==6 && p.termEquals(id, 0) && p.parseInt(1)>=0 && p.termEquals("PREPARED_NOT_SEARCHED", 5),
					"Incomplete or foreign bundled preparation summary");
				count=p.parseLong(2);
				bases=p.parseLong(3);
				records=p.parseLong(4);
			}else{
				binding=oneRow(new File(prepared, "binding.tsv").getPath(), BINDING_HEADER);
				p.set(binding);
				require(p.terms()==11 && p.termEquals(id, 0) && p.termEquals("PREPARATION_BOUND", 10), "Incomplete or foreign preparation binding");
				fasta=text(p, 1);
				fastaSha80=p.parseString(2);
				mapping=text(p, 3);
				mapSha80=p.parseString(4);
				members=text(p, 5);
				membersSha80=p.parseString(6);
				count=p.parseLong(7);
				bases=p.parseLong(8);
				records=p.parseLong(9);
				require(fasta.equals(new File(prepared, "batch.fa").getPath()) && mapping.equals(new File(prepared, "targets.tsv").getPath())
					&& members.equals(new File(prepared, "members.tsv").getPath()), "Preparation binding must name this batch's actual artifacts");
				require(SeedOffsetTable.validSha80(fastaSha80) && SeedOffsetTable.validSha80(mapSha80) && SeedOffsetTable.validSha80(membersSha80)
					&& count>0 && count<=catalog.size() && bases>0 && records>0, "Preparation needs declared sha80 identities and positive dimensions");
				p.set(oneRow(new File(prepared, "summary.tsv").getPath(), PREPARED_HEADER));
				require(p.terms()==6 && p.termEquals(id, 0) && p.parseInt(1)>=0 && p.parseLong(2)==count
					&& p.parseLong(3)==bases && p.parseLong(4)==records && p.termEquals("PREPARED_NOT_SEARCHED", 5), "Preparation summary differs from binding");
			}
			require(count>0 && count<=catalog.size() && bases>0 && records>0, "Prepared batch dimensions must be positive and bounded by the catalog");
			final LongIntMap seen=new LongIntMap(16);
			long memberBases=0, memberRecords=0;
			try(Input in=new Input(members)){
				RrnaResourceIO.header(in.reader, NcrnaEvidenceCatalog.CENSUS_HEADER, members);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					p.set(row);
					require(p.terms()==5, "Prepared membership needs the exact census schema");
					final Genome genome=catalog.get(p.parseString(0));
					require(genome!=null && text(p, 1).equals(genome.fasta) && p.parseLong(2)==genome.bytes
						&& p.parseLong(3)==genome.bases && p.parseLong(4)==genome.records, "Prepared member differs from the exact source census");
					final long number=Long.parseLong(genome.id);
					require(number>0 && !seen.contains(number), "Repeated input in prepared batch");
					seen.put(number, genomes.size());
					genomes.add(genome);
					memberBases=Math.addExact(memberBases, genome.bases);
					memberRecords=Math.addExact(memberRecords, genome.records);
				}
			}
			require(genomes.size()==count && memberBases==bases && memberRecords==records, "Prepared members do not conserve the batch population");
		}
		final String id, key, bundle, prepared, results, fasta, mapping, members, fastaSha80, mapSha80, membersSha80;
		/** Read location only; never substituted into original command bindings or evidence keys. */
		final String physicalResults;
		final byte[] binding;
		final long expectedHits, bases, records;
		final ArrayList<Genome> genomes=new ArrayList<Genome>();
	}

	/** Actual bundled searches retain three provenance rows, not an old-format binding. */
	private static String[] readInputReceipt(final String path, final String prepared){
		final String[] names={"batch.fa", "targets.tsv", "members.tsv"};
		final String[] hashes=new String[names.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "role\tpath\tsha80", path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==3, "Native input receipt requires role, path and sha80");
				int slot=-1;
				for(int i=0; i<names.length; i++){if(p.termEquals(names[i], 0)){slot=i;}}
				require(slot>=0 && hashes[slot]==null && p.termEquals(new File(prepared, names[slot]).getPath(), 1),
					"Unknown, repeated or foreign native input receipt row");
				hashes[slot]=p.parseString(2);
				require(SeedOffsetTable.validSha80(hashes[slot]), "Native input receipt must retain the declared sha80");
			}
		}
		for(final String hash:hashes){require(hash!=null, "Missing native input receipt row");}
		return hashes;
	}

	private static final class Mapping{
		Mapping(final int input_, final String synthetic_, final String original_){
			input=input_;
			synthetic=synthetic_;
			original=original_;
		}
		final int input;
		final String synthetic, original;
	}

	/** Reused token offsets for the few whitespace-separated native accounting lines. */
	private static final class Words{
		void set(final byte[] row){
			data=row;
			count=0;
			int pos=0;
			while(pos<row.length && count<start.length){
				while(pos<row.length && row[pos]<=32){pos++;}
				if(pos==row.length){break;}
				start[count]=pos;
				while(pos<row.length && row[pos]>32){pos++;}
				end[count++]=pos;
			}
		}
		boolean eq(final int field, final String text){
			if(field>=count || end[field]-start[field]!=text.length()){return false;}
			for(int i=0; i<text.length(); i++){if(data[start[field]+i]!=text.charAt(i)){return false;}}
			return true;
		}
		long number(final int field, final int left, final int right){
			require(field<count && end[field]-start[field]>left+right, "Missing native numeric accounting field");
			return Parse.parseLong(data, start[field]+left, end[field]-right);
		}
		long clen(final int field){
			require(field<count && end[field]-start[field]>7, "Missing native query CLEN");
			final String prefix="[CLEN=";
			for(int i=0; i<prefix.length(); i++){require(data[start[field]+i]==prefix.charAt(i), "Malformed native query CLEN");}
			require(data[end[field]-1]==']', "Malformed native query CLEN terminator");
			return number(field, 6, 1);
		}
		String text(final int field){
			require(field<count, "Missing native metadata token");
			return new String(data, start[field], end[field]-start[field], StandardCharsets.UTF_8);
		}
		byte[] data;
		int count;
		final int[] start=new int[8], end=new int[8];
	}

	private static boolean same(final LineParser1 p, final int field, final ByteBuilder expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length()){return false;}
		for(int i=0; i<expected.length(); i++){if(p.line()[p.a()+i]!=expected.array[i]){return false;}}
		return true;
	}
	private static byte[] oneRow(final String path, final String header){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, header, path);
			final byte[] row=in.reader.nextLine();
			require(row!=null && in.reader.nextLine()==null, "Expected one complete metadata row: "+path);
			return row;
		}
	}
	private static String line(final String path){
		try(Input in=new Input(path)){
			final byte[] row=in.reader.nextLine();
			require(row!=null && in.reader.nextLine()==null, "Expected one exit record: "+path);
			return new String(row, StandardCharsets.UTF_8);
		}
	}
	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private final Task task;
	private final Roster roster;
	private final double searchSpaceMb;
	private final Evidence[] evidence;
	private final long[] expected;
	private final long[] observed;
	private final double[] reportBits;
	private final LongIntMap inputIndex=new LongIntMap(16);
	private final ObjectMap<Target,Mapping> targets=new ObjectMap<Target,Mapping>(Target.class, Mapping.class);
	private final Target probe=new Target();
	static final String BINDING_HEADER="batch_id\tfasta\tfasta_sha80\tmapping\tmapping_sha80\tmembers\tmembers_sha80\tgenomes\tnucleotide_bases\tfasta_records\tstatus";
	static final String COMPLETION_HEADER="batch_id\tmodel\tFASTA_records\tinput_bases\tsearched_residues\treported_hits\toriginal_GA_bits\treport_bits\tstatus";
	static final String BUNDLED_COMPLETION_HEADER="batch_id\tmodel\tFASTA_records\tinput_bases\tsearched_residues\treported_hits\tstatus";
	static final String PREPARED_HEADER="batch_id\ttask_id\tgenomes\tnucleotide_bases\tfasta_records\tstatus";
}
