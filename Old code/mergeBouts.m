function Bouts = mergeBouts(Bouts, mergeFrameSep )
    b = 1;
    while b < numel(Bouts)%for b = 1:numel(bout)-1
        if Bouts(b+1).frames(1) - Bouts(b).frames(end) <= mergeFrameSep
            Bouts(b).frames = Bouts(b).frames(1):Bouts(b+1).frames(end);
            Bouts(b).Nframes = numel( Bouts(b).frames );
            Bouts(b+1) = [];
        end
        b = b + 1;
    end
end