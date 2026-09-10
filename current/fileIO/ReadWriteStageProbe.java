package fileIO;

/** Linux diagnostic for finishWriting and buffered-writer error propagation.
 * Prints normal/close/buffered error flags; fixed behavior is false/true/true.
 * Intentionally writes a short record to /dev/full; no ordinary files are created. */
import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import java.io.*;
public final class ReadWriteStageProbe {
 public static void main(String[] args) throws Exception {
  if(!new File("/dev/full").exists())throw new IllegalStateException("This diagnostic requires Linux /dev/full");
  ByteArrayOutputStream good=new ByteArrayOutputStream();
  good.write("valid".getBytes("UTF-8"));
  boolean normal=ReadWrite.finishWriting(null,good,null,false);
  if(normal||!"valid".equals(good.toString("UTF-8")))throw new AssertionError("normal output changed");
  OutputStream broken=new OutputStream(){public void write(int b){} public void close() throws IOException{throw new IOException("intentional close failure");}};
  boolean closeError=ReadWrite.finishWriting(null,broken,null,false);
  ByteStreamWriter writer=new ByteStreamWriter("/dev/full",true,false,true);
  writer.start();writer.print("intentional buffered write failure\n");
  boolean writerError=writer.poisonAndWait();
  System.out.println("normal_error="+normal+" close_error="+closeError+" buffered_writer_error="+writerError);
 }
}
