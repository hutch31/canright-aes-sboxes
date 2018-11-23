module compare_output;

  reg [8:0] sbin;
  wire [7:0] canright_out;
  wire [7:0] reference_out;

  initial
    begin
      $dumpfile("sbox.vcd");
      $dumpvars;
      for (sbin=0; sbin < 256; sbin=sbin+1)
        begin
          #1;
          if (canright_out !== reference_out)
            $display("Miscompare canright=%x reference=%x", canright_out, reference_out);
        end
      $finish;
    end

  canright_sbox can0 (.A(sbin[7:0]), .encrypt(1'b1), .Q(canright_out));
  reference_sbox ref0 (.index(sbin[7:0]), .sbox(reference_out));

endmodule
