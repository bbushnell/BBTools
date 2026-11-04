package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec.Model;
import shared.Tools;
import structures.ByteBuilder;

/** Fixed ordered consensus records with optional aligned HBM records.
 * Arrays are owned here. Individual records are immutable after loading, so a
 * snapshot copies only the outer arrays and rollback restores whole records.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaModelPairs{

	/** Selects a consensus subset; supplied HBMs must match the entire original axis. */
	NcrnaModelPairs(final Library source, final String hbm, final boolean[] selected){
		assert(source.names.length==selected.length) : "Paired selection uses the original consensus index axis";
		int count=0;
		for(final boolean keep:selected){if(keep){count++;}}
		if(count==0){throw new IllegalArgumentException("A library subset needs at least one selected model");}
		names=new String[count];
		headers=new String[count];
		models=hbm==null ? null : new Model[count];
		unpairedBases=hbm==null ? new byte[count][] : null;
		try(NucleotideHbmCodec.Reader reader=hbm==null ? null : new NucleotideHbmCodec.Reader(hbm)){
			int next=0;
			for(int i=0; i<selected.length; i++){
				final Model model=reader==null ? null : reader.next();
				if(reader!=null){
					if(model==null){throw new IllegalArgumentException("Fewer HBMs than input consensuses");}
					Tools.toUpperCase(model.bases);
					if(!name(model).equals(source.names[i]) || !Arrays.equals(model.bases, source.bases[i])){
						throw new IllegalArgumentException("Consensus/HBM pair differs at model "+source.names[i]);
					}
				}
				if(selected[i]){
					names[next]=source.names[i];
					headers[next]=source.headers[i];
					if(models==null){
						unpairedBases[next]=source.bases[i];
					}else{
						models[next]=model;
					}
					next++;
				}
			}
			if(reader!=null && reader.next()!=null){throw new IllegalArgumentException("More HBMs than input consensuses");}
		}
	}

	private NcrnaModelPairs(final NcrnaModelPairs source){
		names=source.names;
		headers=source.headers.clone();
		models=source.models==null ? null : source.models.clone();
		unpairedBases=source.unpairedBases;
	}

	/** Snapshots the mutable selection without duplicating retained sequence/count data. */
	NcrnaModelPairs copy(){return new NcrnaModelPairs(this);}

	/** Creates a borrowed sequence view; no pair may change until its scan returns. */
	Library library(){
		if(models==null){return new Library(names, unpairedBases);}
		final byte[][] bases=new byte[models.length][];
		for(int i=0; i<bases.length; i++){bases[i]=models[i].bases;}
		return new Library(names, bases);
	}

	/** Restores both the full FASTA header and every HBM field at one stable index. */
	boolean restore(final NcrnaModelPairs baseline, final int model){
		assert(Arrays.equals(names, baseline.names)) : "Rollback requires unchanged ordered model identities";
		if(same(baseline, model)){return false;}
		headers[model]=baseline.headers[model];
		models[model]=baseline.models[model];
		return true;
	}

	/** Compares the complete paired record, including header metadata and positional counts. */
	boolean same(final NcrnaModelPairs other, final int i){
		assert(models!=null && other.models!=null) : "Regeneration rollback compares complete consensus/HBM pairs";
		return headers[i].equals(other.headers[i]) && Arrays.equals(models[i].header, other.models[i].header)
			&& Arrays.equals(models[i].bases, other.models[i].bases) && Arrays.equals(models[i].counts, other.models[i].counts);
	}

	/** Writes selected pairs once, preserving requested HBM encoding and compression. */
	void write(final File directory, final String modelFilename){
		assert((models==null)==(modelFilename==null)) : "HBM export requires paired records and an explicit codec filename together";
		final ByteBuilder row=new ByteBuilder();
		try(Output fasta=new Output(new File(directory, "consensus.fa").getPath());
				NucleotideHbmCodec.Writer hbm=models==null ? null : new NucleotideHbmCodec.Writer(new File(directory, modelFilename).getPath(),
					NucleotideHbmCodec.isCompact(modelFilename), false)){
			for(int i=0; i<names.length; i++){
				row.clear();
				fasta.writer.print(row.append('>').append(headers[i]).nl().append(models==null ? unpairedBases[i] : models[i].bases).nl());
				if(hbm!=null){hbm.add(models[i]);}
			}
			if(hbm!=null){hbm.finish();}
		}
	}

	/** HBM identifiers precede optional tab-separated count metadata. */
	static String name(final Model model){return NcrnaSupport.recordName(header(model));}

	/** FASTA metadata ends where the HBM's tab-separated count fields begin. */
	static String header(final Model model){
		return new String(model.header, 1, headerEnd(model)-1, StandardCharsets.US_ASCII);
	}

	/** Renames only the identifier, retaining description and all HBM metadata bytes. */
	static Model rename(final Model model, final String name){
		final int end=headerEnd(model);
		int nameEnd=1;
		while(nameEnd<end && !Character.isWhitespace((char)model.header[nameEnd])){nameEnd++;}
		final byte[] header=new ByteBuilder().append('>').append(name)
			.append(model.header, nameEnd, model.header.length-nameEnd).toBytes();
		return new Model(header, model.bases, model.counts);
	}

	/** Shared delimiter scan for names, FASTA headers and paired-record renaming. */
	private static int headerEnd(final Model model){
		assert(model.header.length>1 && model.header[0]=='>') : "NucleotideHbmCodec returns a nonempty FASTA-style model header";
		int end=1;
		while(end<model.header.length && model.header[end]!='\t'){end++;}
		return end;
	}

	final String[] names, headers;
	final Model[] models;
	/** Borrowed immutable sequences for a consensus-only subset; null when paired. */
	private final byte[][] unpairedBases;
}
