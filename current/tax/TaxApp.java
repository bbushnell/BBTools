package tax;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 * Simple Swing GUI application for taxonomy-related operations.
 * Creates a basic window with an interactive button interface.
 * @author Brian Bushnell
 */
public class TaxApp {
	
	public static void main(final String[] args) {
        final JFrame parent = new JFrame();
        JButton button = new JButton();

        button.setText("Button text");
        parent.add(button);
        parent.pack();
        parent.setVisible(true);

        button.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                String name = JOptionPane.showInputDialog(parent,
                        "Prompt", null);
            }
        });
    }
	
}