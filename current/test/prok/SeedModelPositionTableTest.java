package test.prok;

import java.io.File;
import java.util.Arrays;

import aligner.CovarianceModel;
import aligner.CovarianceModelConsensusMap;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.ObjectMap;
import parse.LineParser1;
import prok.RrnaResourceIO;
import prok.SeedModelPositionResource;
import prok.SeedModelPositionTable;
import structures.ByteBuilder;

/** Hand-constructed global maps test exact annotation without scoring/alignment.
 * @author Brian Bushnell, Raiden */
public final class SeedModelPositionTableTest {

	public static void main(final String[] args){
		check(args.length==1, "Pass one new fixture directory");
		final File dir=new File(args[0]);
		check(dir.mkdir(), "Fixture directory must be new and writable: "+dir);
		final CovarianceModel cm=pairModel();
		final byte[] a=new byte[37], g=new byte[17], rfA=new byte[37], rfG=new byte[17];
		Arrays.fill(a, (byte)'A');Arrays.fill(g, (byte)'G');
		Arrays.fill(rfA, (byte)'.');Arrays.fill(rfG, (byte)'.');
		rfA[17]='x';rfA[19]='x';rfG[0]='x';rfG[16]='x';
		final byte[][] bases={a, g};
		final CovarianceModelConsensusMap[] maps={
			new CovarianceModelConsensusMap(cm, a, a, rfA), new CovarianceModelConsensusMap(cm, g, g, rfG)};
		long gKey=0;
		for(int i=0; i<17; i++){gKey=(gKey<<2)|2;}
		final LongHashSet seeds=new LongHashSet(8);
		seeds.add(0);seeds.add(gKey);seeds.add(SeedModelPositionTable.MASK);
		final String[] names={"repeat_A", "single_G"};
		final String[] hashes={"00000000000000000001", "00000000000000000002",
			"00000000000000000003", "00000000000000000004"};
		final int[] lengths={37, 17}, sites={1};
		final SeedModelPositionTable table=SeedModelPositionTable.build(seeds, names, bases, maps, sites, hashes);
		checkPositions(table, gKey);
		roundTrip(dir, table, seeds, names, lengths, sites, hashes);
		subset(dir, table, seeds, names, bases, gKey);
		System.out.println("SEED_MODEL_POSITION_TEST_PASS repeat unmapped orientation interval_sides "
			+"terminal_insertion provenance conservation subset alignments=0 cmScores=0");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Checks            ----------------*/
	/*--------------------------------------------------------------*/

	private static void checkPositions(final SeedModelPositionTable table, final long gKey){
		check(table.seedCount()==3 && table.get(0).size()==21 && table.get(0).ambiguous(),
			"All overlapping exact positions must survive, rather than keeping one convenient anchor");
		check(table.get(gKey).size()==1 && !table.get(gKey).ambiguous()
			&& table.get(SeedModelPositionTable.MASK).size()==0,
			"Forward-only annotation retains absent T17 as unmapped, not as reverse-complement support from A17");
		check(table.side(table.get(0).get(1), 0)==-1 && table.side(table.get(0).get(2), 0)==0
			&& table.side(table.get(0).get(19), 0)==1,
			"The RF boundary maps to query interval18..19; only whole seeds beyond that interval support a side");
		check(table.get(0).get(0).rfStartLow==0 && table.get(0).get(0).rfEndHigh==0,
			"A seed entirely inside a terminal insertion maps to that boundary, never fabricated model residues");
	}

	private static void roundTrip(final File dir, final SeedModelPositionTable table, final LongHashSet seeds,
			final String[] names, final int[] lengths, final int[] sites, final String[] hashes){
		final String path=new File(dir, "positions.tsv").toString();
		table.write(path);
		final SeedModelPositionTable loaded=SeedModelPositionTable.load(path, seeds, names, lengths, 2, sites, hashes);
		check(loaded.get(0).size()==21 && loaded.get(SeedModelPositionTable.MASK).size()==0
			&& loaded.side(loaded.get(0).get(19), 0)==1,
			"The native loader preserves repeated positions, the unmapped sentinel and projected-site semantics");
		final String[] wrong=hashes.clone();
		wrong[1]=SUBSET_HASH;
		boolean rejected=false;
		try{SeedModelPositionTable.load(path, seeds, names, lengths, 2, sites, wrong);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A resource from another consensus generation must fail before it can place seeds");
		final String missing=new File(dir, "missing.tsv").toString();
		copyWithoutLastLine(path, missing);
		rejected=false;
		try{SeedModelPositionTable.load(missing, seeds, names, lengths, 2, sites, hashes);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Unmapped keys are still required resource members, including the last key");
	}

	private static void subset(final File dir, final SeedModelPositionTable source, final LongHashSet seeds,
			final String[] names, final byte[][] bases, final long gKey){
		final SeedModelPositionTable onlyG=source.subset(new int[]{1}, SUBSET_HASH);
		check(onlyG.seedCount()==3 && onlyG.modelCount()==1 && onlyG.modelName(0).equals("single_G")
			&& onlyG.get(0).size()==0 && onlyG.get(gKey).size()==1 && onlyG.get(gKey).get(0).model==0,
			"Removing the first model remaps the second to zero while its removed-only seed remains explicitly unmapped");
		final SeedModelPositionTable onlyA=source.subset(new int[]{0}, SUBSET_HASH);
		check(onlyA.get(0).size()==21 && onlyA.get(0).ambiguous() && onlyA.get(gKey).size()==0
			&& onlyA.side(onlyA.get(0).get(1), 0)==-1 && onlyA.side(onlyA.get(0).get(2), 0)==0
			&& onlyA.side(onlyA.get(0).get(19), 0)==1 && source.modelCount()==2 && source.get(gKey).get(0).model==1,
			"Subsetting preserves all repeat postings and uncertain site intervals without mutating the source");
		for(int[] bad:new int[][]{new int[0], {1, 0}, {0, 0}, {2}, {-1}}){
			boolean rejected=false;
			try{source.subset(bad, SUBSET_HASH);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Invalid subset indices must fail before constructing an incompatible model map: "+Arrays.toString(bad));
		}
		boolean rejected=false;
		try{source.subset(new int[]{0}, "invalid");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "New consensus provenance must satisfy the table's sha80 contract");
		subsetFiles(dir, source, seeds, names, bases, gKey);
	}

	private static void subsetFiles(final File dir, final SeedModelPositionTable source, final LongHashSet seeds,
			final String[] names, final byte[][] bases, final long gKey){
		final String envelope=new File(dir, "source.resource.tsv").toString();
		final String modelFile=new File(dir, "models.fa").toString();
		final String selectedFile=new File(dir, "selected.fa").toString();
		final String seedFile=new File(dir, "seeds.fa").toString();
		SeedModelPositionResource.write(envelope, source, bases);
		final ByteBuilder text=new ByteBuilder();
		writeFixture(modelFile, text.append('>').append(names[0]).nl().append(bases[0]).nl()
			.append('>').append(names[1]).nl().append(bases[1]).nl());
		writeFixture(selectedFile, text.clear().append(">single_G\n").append(bases[1]).nl());
		writeFixture(seedFile, text.clear().append(">A\nAAAAAAAAAAAAAAAAA\n>G\nGGGGGGGGGGGGGGGGG\n>T\nTTTTTTTTTTTTTTTTT\n"));
		final String out=new File(dir, "selected.resource.tsv").toString();
		SeedModelPositionResource.main(new String[]{"positions="+envelope, "models="+modelFile, "seeds="+seedFile,
			"subsetmodels="+selectedFile, "subsetconsensussha80="+SUBSET_HASH, "out="+out});
		final SeedModelPositionTable loaded=SeedModelPositionResource.load(out, seeds, new String[]{"single_G"}, new byte[][]{bases[1]});
		check(loaded.modelCount()==1 && loaded.get(gKey).get(0).model==0 && loaded.get(0).size()==0
			&& loaded.get(SeedModelPositionTable.MASK).size()==0,
			"The actual subset CLI must round-trip a remapped model and both preexisting and newly unmapped keys");
		checkBindings(out);
		final String changed=new File(dir, "changed.fa").toString();
		writeFixture(changed, text.clear().append(">single_G\nCGGGGGGGGGGGGGGGG\n"));
		final File forbidden=new File(dir, "must_not_exist.tsv");
		boolean rejected=false;
		try{
			SeedModelPositionResource.main(new String[]{"positions="+envelope, "models="+modelFile, "seeds="+seedFile,
				"subsetmodels="+changed, "subsetconsensussha80="+SUBSET_HASH, "out="+forbidden});
		}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !forbidden.exists(),
			"A changed consensus is not a subset; reject it before writing inherited positions against different bases");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Explicit S -> MP -> E fixture; no reflective dependency on another test. */
	private static CovarianceModel pairModel(){
		final ObjectMap<String, String> headers=new ObjectMap<String, String>(String.class, String.class);
		headers.put("QDBBETA1", "1e-7");headers.put("QDBBETA2", "1e-15");
		final CovarianceModel model=new CovarianceModel("pair", "TEST", 2, 20, 3, 3, headers,
			new float[]{.25f, .25f, .25f, .25f});
		model.type[0]=CovarianceModel.S;model.type[1]=CovarianceModel.MP;model.type[2]=CovarianceModel.E;
		for(int i=0; i<3; i++){
			model.node[i]=i;model.childFirst[i]=(i<2 ? i+1 : -1);model.childCount[i]=(i<2 ? 1 : 0);
			model.transitionScore[i]=(i<2 ? new float[]{0} : new float[0]);
			model.qdb[i]=(i<2 ? new int[]{0, 0, 20, 20} : new int[]{0, 0, 0, 0});
			model.emissionScore[i]=new float[i==1 ? 16 : 0];
			Arrays.fill(model.emissionScore[i], -1);
		}
		model.consensusLeft[1]=1;model.consensusRight[1]=2;model.emissionScore[1][3]=8;
		return model;
	}

	private static void writeFixture(final String path, final ByteBuilder text){
		check(text.length()>0, "A fixture must contain its declared rows before the writer starts");
		final ByteStreamWriter out=RrnaResourceIO.writer(path);
		try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}

	private static void copyWithoutLastLine(final String source, final String target){
		final ByteFile in=RrnaResourceIO.open(source);
		final ByteStreamWriter out=RrnaResourceIO.writer(target);
		try{
			byte[] previous=in.nextLine();
			check(previous!=null, "The valid position table must contain a header before removing its last row");
			for(byte[] next; (next=in.nextLine())!=null;){out.println(previous);previous=next;}
		}finally{
			boolean error=in.close();
			error|=out.poisonAndWait();
			check(!error, "Both fixture streams must close successfully before testing missing-key rejection");
		}
	}

	private static void checkBindings(final String path){
		final String[] hashes={"00000000000000000001", SUBSET_HASH, "00000000000000000003", "00000000000000000004"};
		final ByteFile in=RrnaResourceIO.open(path);
		final LineParser1 parser=new LineParser1('\t');
		int found=0;
		try{
			for(byte[] line; (line=in.nextLine())!=null;){
				parser.set(line);
				if(!parser.termEquals("#Binding", 0)){continue;}
				check(parser.terms()==3 && found<hashes.length
					&& parser.termEquals(SeedModelPositionTable.BINDINGS[found], 1) && parser.termEquals(hashes[found], 2),
					"Only consensus provenance changes; ordered seed, anchor and CM bindings must survive subsetting");
				found++;
			}
		}finally{check(!in.close(), "The serialized binding fixture must be read completely");}
		check(found==hashes.length, "Every original provenance binding must remain in the subset resource");
	}

	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String SUBSET_HASH="00000000000000000005";
}
