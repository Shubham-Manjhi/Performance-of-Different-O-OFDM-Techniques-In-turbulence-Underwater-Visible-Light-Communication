function [out,in] = vlc_led_filter(in, led_filter, led)
    out = led_filter(in,led);
end